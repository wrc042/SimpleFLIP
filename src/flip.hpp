#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <thread>
#include <vector>

class FLIPSolver {
  public:
    FLIPSolver(int resolution, double time_interval, double num_marker)
        : _resolution{resolution}, _time_interval{time_interval},
          _num_marker(num_marker) {
        std::random_device rd;
        auto seed{rd()};
        _rand_engine = std::default_random_engine(seed);

        _num_particle = 0;
        _grid_spacing = 1.0 / _resolution;

        _velocityu.resize(_resolution * (_resolution + 1), 0);
        _velocityu_.resize(_resolution * (_resolution + 1), 0);
        _weightu.resize(_resolution * (_resolution + 1), 0);
        _velocityv.resize(_resolution * (_resolution + 1), 0);
        _velocityv_.resize(_resolution * (_resolution + 1), 0);
        _weightv.resize(_resolution * (_resolution + 1), 0);
        _marker.resize(_resolution * _resolution);
        _density.resize(_resolution * _resolution);

        for (int i = 0; i < _resolution * _resolution; i++) {
            _marker[i] = AIR;
        }
        for (int i = 0; i < _resolution; i++) {
            _marker[i * _resolution + 0] = SOLID;
            _marker[i * _resolution + 1] = SOLID;
            _marker[i * _resolution + _resolution - 1] = SOLID;
            _marker[i * _resolution + _resolution - 2] = SOLID;
            _marker[(0) * _resolution + i] = SOLID;
            _marker[(1) * _resolution + i] = SOLID;
            _marker[(_resolution - 1) * _resolution + i] = SOLID;
            _marker[(_resolution - 2) * _resolution + i] = SOLID;
        }
    }

    double radius() { return _grid_spacing / sqrt(_num_marker); }
    void set_reset_buffer(
        const std::function<void(const std::vector<Eigen::Vector2d> &)> &func) {
        _reset_buffer = func;
    }
    void set_update_buffer(
        const std::function<void(const std::vector<Eigen::Vector2d> &)> &func) {
        _update_buffer = func;
    }
    void step() {
        clear_grid();
        particle2grid();
        gravity(_time_interval / _substep);
        projection();
        grid2particle();
        advection();
    }
    double get_radius() { return _grid_spacing / sqrt(_num_marker); }
    void run(const bool &is_closed) {
        if (_reset_buffer != NULL) {
            _reset_buffer(_particle_position);
        }
        while (!is_closed) {
            double st = omp_get_wtime();
            for (int i = 0; i < _substep; i++) {
                step();
            }
            if (_update_buffer != NULL) {
                _update_buffer(_particle_position);
            }
            double et = omp_get_wtime();
            if ((1.0 / _max_fps - (et - st)) > 0) {
                std::this_thread::sleep_for(std::chrono::milliseconds(
                    (int)(1000 * (1.0 / _max_fps - (et - st)))));
            }
        }
    }
    void add_particle(const std::function<double(Eigen::Vector2d)> &sdf) {
        std::uniform_real_distribution<double> dist(0.0, _grid_spacing);
        for (int i = 0; i < _resolution; i++) {
            for (int j = 0; j < _resolution; j++) {
                Eigen::Vector2d cell_corner(i * _grid_spacing,
                                            j * _grid_spacing);
                for (int k = 0; k < _num_marker; k++) {
                    double offsetx = dist(_rand_engine);
                    double offsety = dist(_rand_engine);
                    Eigen::Vector2d pos =
                        cell_corner + Eigen::Vector2d(offsetx, offsety);
                    if (sdf(pos) < 0.0) {
                        _num_particle += 1;
                        _particle_position.push_back(pos);
                        _particle_velocity.push_back(Eigen::Vector2d::Zero());
                    }
                }
            }
        }
    }

  private:
    void clear_grid() {
#pragma omp parallel for
        for (int i = 0; i < _resolution * (_resolution + 1); i++) {
            _velocityu[i] = 0.0;
            _weightu[i] = 0.0;
            _velocityv[i] = 0.0;
            _weightv[i] = 0.0;
        }
    }

    void particle2grid() {
        double inv_grid_spacing = 1.0 / _grid_spacing;
        auto get_idx = [inv_grid_spacing](double x, double y) {
            int xidx = floor(x * inv_grid_spacing);
            int yidx = floor(y * inv_grid_spacing);
            return std::array<int, 2>{xidx, yidx};
        };
        auto get_frac = [inv_grid_spacing](double x, double y, int xidx,
                                           int yidx) {
            double fracx = x * inv_grid_spacing - xidx;
            double fracy = y * inv_grid_spacing - yidx;
            return std::array<double, 4>{fracx * fracy, (1 - fracx) * fracy,
                                         fracx * (1 - fracy),
                                         (1 - fracx) * (1 - fracy)};
        };
        for (int i = 0; i < _num_particle; i++) {
            std::array<int, 4> offsetx = {0, 1, 0, 1};
            std::array<int, 4> offsety = {0, 0, 1, 1};

            auto idxu =
                get_idx(_particle_position[i].x(),
                        _particle_position[i].y() - 0.5 * _grid_spacing);
            auto fracu =
                get_frac(_particle_position[i].x(),
                         _particle_position[i].y() - 0.5 * _grid_spacing,
                         idxu[0], idxu[1]);
            auto idxv = get_idx(_particle_position[i].x() - 0.5 * _grid_spacing,
                                _particle_position[i].y());
            auto fracv =
                get_frac(_particle_position[i].x() - 0.5 * _grid_spacing,
                         _particle_position[i].y(), idxv[0], idxv[1]);

            for (int j = 0; j < 4; j++) {
                int tmpidx = 0;
                tmpidx = (idxu[0] + offsetx[j]) * _resolution +
                         (idxu[1] + offsety[j]);
                _velocityu[tmpidx] += _particle_velocity[i].x() * fracu[j];
                _weightu[tmpidx] += fracu[j];

                tmpidx = (idxv[0] + offsetx[j]) * (_resolution + 1) +
                         (idxv[1] + offsety[j]);

                _velocityv[tmpidx] += _particle_velocity[i].y() * fracv[j];
                _weightv[tmpidx] += fracv[j];
            }
        }
#pragma omp parallel for
        for (int i = 0; i < _resolution * (_resolution + 1); i++) {
            if (_weightu[i] > 1e-6) {
                _velocityu[i] /= _weightu[i];
            }
            if (_weightv[i] > 1e-6) {
                _velocityv[i] /= _weightv[i];
            }
            _velocityu_[i] = _velocityu[i];
            _velocityv_[i] = _velocityv[i];
        }
    }
    void grid2particle() {
        double inv_grid_spacing = 1.0 / _grid_spacing;
        auto get_idx = [inv_grid_spacing](double x, double y) {
            int xidx = floor(x * inv_grid_spacing);
            int yidx = floor(y * inv_grid_spacing);
            return std::array<int, 2>{xidx, yidx};
        };
        auto get_frac = [inv_grid_spacing](double x, double y, int xidx,
                                           int yidx) {
            double fracx = x * inv_grid_spacing - xidx;
            double fracy = y * inv_grid_spacing - yidx;
            return std::array<double, 4>{fracx * fracy, (1 - fracx) * fracy,
                                         fracx * (1 - fracy),
                                         (1 - fracx) * (1 - fracy)};
        };
#pragma omp parallel for
        for (int i = 0; i < _num_particle; i++) {
            std::array<int, 4> offsetx = {0, 1, 0, 1};
            std::array<int, 4> offsety = {0, 0, 1, 1};

            double vu = 0.0;
            double vv = 0.0;
            double dvu = 0.0;
            double dvv = 0.0;

            auto idxu =
                get_idx(_particle_position[i].x(),
                        _particle_position[i].y() - 0.5 * _grid_spacing);
            auto fracu =
                get_frac(_particle_position[i].x(),
                         _particle_position[i].y() - 0.5 * _grid_spacing,
                         idxu[0], idxu[1]);
            auto idxv = get_idx(_particle_position[i].x() - 0.5 * _grid_spacing,
                                _particle_position[i].y());
            auto fracv =
                get_frac(_particle_position[i].x() - 0.5 * _grid_spacing,
                         _particle_position[i].y(), idxv[0], idxv[1]);

            for (int j = 0; j < 4; j++) {
                int tmpidx = 0;
                tmpidx = (idxu[0] + offsetx[j]) * _resolution +
                         (idxu[1] + offsety[j]);
                vu += fracu[j] * _velocityu[tmpidx];
                dvu += fracu[j] * (_velocityu[tmpidx] - _velocityu_[tmpidx]);

                tmpidx = (idxv[0] + offsetx[j]) * (_resolution + 1) +
                         (idxv[1] + offsety[j]);
                vv += fracv[j] * _velocityv[tmpidx];
                dvv += fracv[j] * (_velocityv[tmpidx] - _velocityv_[tmpidx]);
            }
            double vflipu = _particle_velocity[i].x() + dvu;
            double vflipv = _particle_velocity[i].y() + dvv;
            _particle_velocity[i].x() =
                _flip_weight * vflipu + (1 - _flip_weight) * vu;
            _particle_velocity[i].y() =
                _flip_weight * vflipv + (1 - _flip_weight) * vv;
        }
    }

    void gravity(double time_interval) {
#pragma omp parallel for
        for (int i = 0; i < _resolution * (_resolution + 1); i++) {
            _velocityv[i] -= 9.8 * time_interval;
        }
    }

    void advection() {
#pragma omp parallel for
        for (int i = 0; i < _num_particle; i++) {
            _particle_position[i] += _particle_velocity[i] * _time_interval;
            _particle_position[i].x() =
                (std::max)(_grid_spacing * 2, _particle_position[i].x());
            _particle_position[i].x() =
                (std::min)(1 - _grid_spacing * 2, _particle_position[i].x());
            _particle_position[i].y() =
                (std::max)(_grid_spacing * 2, _particle_position[i].y());
            _particle_position[i].y() =
                (std::min)(1 - _grid_spacing * 2, _particle_position[i].y());
        }
    }

    void projection() {
        double inv_grid_spacing = 1.0 / _grid_spacing;
        auto get_idx = [inv_grid_spacing](double x, double y) {
            int xidx = floor(x * inv_grid_spacing);
            int yidx = floor(y * inv_grid_spacing);
            return std::array<int, 2>{xidx, yidx};
        };
        std::map<std::array<int, 2>, int> grid2mat;
        for (int i = 0; i < _resolution * _resolution; i++) {
            _density[i] = 0.0;
            if (_marker[i] != SOLID) {
                _marker[i] = AIR;
            }
        }
        for (int i = 0; i < _num_particle; i++) {
            auto idx =
                get_idx(_particle_position[i].x(), _particle_position[i].y());
            if (_marker[idx[0] * _resolution + idx[1]] == AIR) {
                _marker[idx[0] * _resolution + idx[1]] = FLUID;
            }
            _density[idx[0] * _resolution + idx[1]] += 1.0;
        }
        int vecsize = 0;
        for (int i = 0; i < _resolution; i++) {
            for (int j = 0; j < _resolution; j++) {
                if (_marker[i * _resolution + j] == FLUID) {
                    grid2mat[{i, j}] = vecsize;
                    vecsize++;
                }
            }
        }
        typedef Eigen::Triplet<double> T;
        std::vector<T> triple_list;
        Eigen::SparseMatrix<double> A(vecsize, vecsize);
        Eigen::VectorXd b(vecsize);
        Eigen::VectorXd p(vecsize);

        int cnt = 0;

        for (int i = 0; i < _resolution; i++) {
            for (int j = 0; j < _resolution; j++) {
                if (_marker[i * _resolution + j] == FLUID) {
                    double btmp =
                        -_grid_spacing *
                        (_velocityu[(i + 1) * _resolution + j] -
                         _velocityu[i * _resolution + j] +
                         _velocityv[i * (_resolution + 1) + j + 1] -
                         _velocityv[i * (_resolution + 1) + j] -
                         _density_factor * _grid_spacing * _grid_spacing *
                             _grid_spacing *
                             (_density[i * _resolution + j] - _num_marker) /
                             _num_marker);
                    double coeff = 0;
                    if (_marker[(i - 1) * _resolution + j] == FLUID) {
                        coeff += 1;
                        triple_list.push_back(T(cnt, grid2mat[{i - 1, j}], -1));
                    } else if (_marker[(i - 1) * _resolution + j] == AIR) {
                        coeff += 1;
                    } else {
                        btmp -= _grid_spacing * _velocityu[i * _resolution + j];
                    }
                    if (_marker[(i + 1) * _resolution + j] == FLUID) {
                        coeff += 1;
                        triple_list.push_back(T(cnt, grid2mat[{i + 1, j}], -1));
                    } else if (_marker[(i + 1) * _resolution + j] == AIR) {
                        coeff += 1;
                    } else {
                        btmp += _grid_spacing *
                                _velocityu[(i + 1) * _resolution + j];
                    }
                    if (_marker[i * _resolution + j - 1] == FLUID) {
                        coeff += 1;
                        triple_list.push_back(T(cnt, grid2mat[{i, j - 1}], -1));
                    } else if (_marker[i * _resolution + j - 1] == AIR) {
                        coeff += 1;
                    } else {
                        btmp -= _grid_spacing *
                                _velocityv[i * (_resolution + 1) + j];
                    }
                    if (_marker[i * _resolution + j + 1] == FLUID) {
                        coeff += 1;
                        triple_list.push_back(T(cnt, grid2mat[{i, j + 1}], -1));
                    } else if (_marker[i * _resolution + j + 1] == AIR) {
                        coeff += 1;
                    } else {
                        btmp += _grid_spacing *
                                _velocityv[i * (_resolution + 1) + j + 1];
                    }
                    triple_list.push_back(T(cnt, grid2mat[{i, j}], coeff));
                    b(cnt) = btmp;
                    cnt++;
                }
            }
        }
        A.setFromTriplets(triple_list.begin(), triple_list.end());
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                 Eigen::Lower | Eigen::Upper>
            cg;
        cg.compute(A);
        p = cg.solve(b);
        if (cg.info() != Eigen::Success) {
            std::cout << "projection solver failed" << std::endl;
            return;
        }

        for (int i = 1; i < _resolution; i++) {
            for (int j = 0; j < _resolution; j++) {
                if (_marker[i * _resolution + j] == FLUID &&
                    _marker[(i - 1) * _resolution + j] == FLUID) {
                    _velocityu[i * _resolution + j] -=
                        (p(grid2mat[{i, j}]) - p(grid2mat[{i - 1, j}])) /
                        _grid_spacing;
                } else if (_marker[i * _resolution + j] == FLUID &&
                           _marker[(i - 1) * _resolution + j] == AIR) {
                    _velocityu[i * _resolution + j] -=
                        p(grid2mat[{i, j}]) / _grid_spacing;
                } else if (_marker[i * _resolution + j] == AIR &&
                           _marker[(i - 1) * _resolution + j] == FLUID) {
                    _velocityu[i * _resolution + j] +=
                        p(grid2mat[{i - 1, j}]) / _grid_spacing;
                }
            }
        }
        for (int i = 0; i < _resolution; i++) {
            for (int j = 1; j < _resolution; j++) {
                if (_marker[i * _resolution + j] == FLUID &&
                    _marker[i * _resolution + j - 1] == FLUID) {
                    _velocityv[i * (_resolution + 1) + j] -=
                        (p(grid2mat[{i, j}]) - p(grid2mat[{i, j - 1}])) /
                        _grid_spacing;
                } else if (_marker[i * _resolution + j] == FLUID &&
                           _marker[i * _resolution + j - 1] == AIR) {
                    _velocityv[i * (_resolution + 1) + j] -=
                        p(grid2mat[{i, j}]) / _grid_spacing;
                } else if (_marker[i * _resolution + j] == AIR &&
                           _marker[i * _resolution + j - 1] == FLUID) {
                    _velocityv[i * (_resolution + 1) + j] +=
                        p(grid2mat[{i, j - 1}]) / _grid_spacing;
                }
            }
        }

        for (int i = 0; i < _resolution; i++) {
            _velocityu[0 * _resolution + i] = 0;
            _velocityu[1 * _resolution + i] = 0;
            _velocityu[(_resolution)*_resolution + i] = 0;
            _velocityu[(_resolution - 1) * _resolution + i] = 0;
            _velocityv[i * (_resolution + 1) + 0] = 0;
            _velocityv[i * (_resolution + 1) + 1] = 0;
            _velocityv[i * (_resolution + 1) + _resolution] = 0;
            _velocityv[i * (_resolution + 1) + _resolution - 1] = 0;
        }
    }

    std::function<void(const std::vector<Eigen::Vector2d> &)> _reset_buffer;
    std::function<void(const std::vector<Eigen::Vector2d> &)> _update_buffer;

    double _time_interval;

    int _resolution;
    double _grid_spacing;
    std::vector<double> _velocityu;
    std::vector<double> _velocityu_;
    std::vector<double> _velocityv;
    std::vector<double> _velocityv_;
    std::vector<double> _weightu;
    std::vector<double> _weightv;
    std::vector<Eigen::Vector2d> _particle_position;
    std::vector<Eigen::Vector2d> _particle_velocity;
    int _num_particle;

    enum Marker { AIR, SOLID, FLUID };
    std::vector<int> _marker;
    std::vector<double> _density;
    // particle num per cell
    int _num_marker;
    std::default_random_engine _rand_engine;

    double _flip_weight = 0.95;
    int _substep = 1;
    double _density_factor = 20;
    double _max_fps = 60;
};
