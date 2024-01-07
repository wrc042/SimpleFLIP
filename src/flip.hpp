#pragma once

#include <Eigen/Dense>
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
        _marker.resize(_resolution * _resolution, AIR);

        for (int i = 0; i < _resolution; i++) {
            _marker[i * _resolution + 0] = SOLID;
            _marker[i * _resolution + 1] = SOLID;
            _marker[i * _resolution + _resolution - 1] = SOLID;
            _marker[i * _resolution + _resolution - 2] = SOLID;
        }
        for (int j = 0; j < _resolution; j++) {
            _marker[(0) * _resolution + j] = SOLID;
            _marker[(1) * _resolution + j] = SOLID;
            _marker[(_resolution - 1) * _resolution + j] = SOLID;
            _marker[(_resolution - 2) * _resolution + j] = SOLID;
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
        gravity(_time_interval);
        grid2particle();
        advection();
    }
    double get_radius() { return _grid_spacing / sqrt(_num_marker); }
    void run(const bool &is_closed) {
        if (_reset_buffer != NULL) {
            _reset_buffer(_particle_position);
        }
        while (!is_closed) {
            step();
            if (_update_buffer != NULL) {
                _update_buffer(_particle_position);
            }
        }
    }
    void add_particle(const std::function<double(Eigen::Vector2d)> &sdf) {
        std::uniform_real_distribution<double> dist(0.0, _grid_spacing);
        for (int i = 0; i < _resolution; i++) {
            for (int j = 0; j < _resolution; j++) {
                Eigen::Vector2d cell_corner(
                    i * _grid_spacing, (_resolution - 1 - j) * _grid_spacing);
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
        int resolution = _resolution;
        int grid_spacing = _grid_spacing;
        auto get_idx = [resolution](double x, double y) {
            int xidx = floor(x * resolution);
            int yidx = floor(y * resolution);
            return std::array<int, 2>{xidx, yidx};
        };
        auto get_frac = [resolution, grid_spacing](double x, double y, int xidx,
                                                   int yidx) {
            double fracx = x * resolution - xidx;
            double fracy = y * resolution - yidx;
            return std::array<double, 4>{fracx * fracy, (1 - fracx) * fracy,
                                         fracx * (1 - fracy),
                                         (1 - fracx) * (1 - fracy)};
        };
        for (int i = 0; i < _num_particle; i++) {
            std::array<int, 4> offsetx = {0, 1, 0, 1};
            std::array<int, 4> offsety = {0, 0, 1, 1};

            auto idxu =
                get_idx(_particle_position[i].x(),
                        _particle_position[i].y() + 0.5 * _grid_spacing);
            auto fracu =
                get_frac(_particle_position[i].x(),
                         _particle_position[i].y() + 0.5 * _grid_spacing,
                         idxu[0], idxu[1]);
            auto idxv = get_idx(_particle_position[i].x() + 0.5 * _grid_spacing,
                                _particle_position[i].y());
            auto fracv =
                get_frac(_particle_position[i].x() + 0.5 * _grid_spacing,
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
            if (_weightu[i] > 0.0) {
                _velocityu[i] /= _weightu[i];
            }
            if (_weightv[i] > 0.0) {
                _velocityv[i] /= _weightv[i];
            }
            _velocityu_[i] = _velocityu[i];
            _velocityv_[i] = _velocityv[i];
        }
    }
    void grid2particle() {
        int resolution = _resolution;
        int grid_spacing = _grid_spacing;
        auto get_idx = [resolution](double x, double y) {
            int xidx = floor(x * resolution);
            int yidx = floor(y * resolution);
            return std::array<int, 2>{xidx, yidx};
        };
        auto get_frac = [resolution, grid_spacing](double x, double y, int xidx,
                                                   int yidx) {
            double fracx = x * resolution - xidx;
            double fracy = y * resolution - yidx;
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
                        _particle_position[i].y() + 0.5 * _grid_spacing);
            auto fracu =
                get_frac(_particle_position[i].x(),
                         _particle_position[i].y() + 0.5 * _grid_spacing,
                         idxu[0], idxu[1]);
            auto idxv = get_idx(_particle_position[i].x() + 0.5 * _grid_spacing,
                                _particle_position[i].y());
            auto fracv =
                get_frac(_particle_position[i].x() + 0.5 * _grid_spacing,
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

    // void projection() {

    //     for (int i = 0; i < _resolution; i++) {
    //         _velocityu[0 * _resolution + i] = 0;
    //         _velocityu[1 * _resolution + i] = 0;
    //         _velocityu[(_resolution)*_resolution + i] = 0;
    //         _velocityu[(_resolution - 1) * _resolution + i] = 0;
    //         _velocityv[i * (_resolution + 1) + 0] = 0;
    //         _velocityv[i * (_resolution + 1) + 1] = 0;
    //         _velocityv[i * (_resolution + 1) + _resolution] = 0;
    //         _velocityv[i * (_resolution + 1) + _resolution - 1] = 0;
    //     }
    // }

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
    // particle num per cell
    int _num_marker;
    std::default_random_engine _rand_engine;

    double _flip_weight = 0.95;
};