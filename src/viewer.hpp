#pragma once

#include "omp.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"

class Viewer {
  public:
    Viewer() { init(); }
    void init() {
        set_window();
        set_scene();
        set_camera();

        polyscope::init();
        set_callback();

        set_ground_after();
        set_scene_after();

        omp_init_lock(&_buflock);
    }
    void show() {
        polyscope::show();
        _is_closed = true;
    }
    const bool &is_closed() { return _is_closed; }
    void set_radius(double radius) { _radius = radius; }

    void reset_buffer(const std::vector<Eigen::Vector2d> &points) {
        std::array<int, 2> bufnum;
        omp_set_lock(&_buflock);
        bufnum = _bufnum;
        omp_unset_lock(&_buflock);

        bufnum[1] = 1;
        bufnum[0] = 0;
        _points[bufnum[1]].clear();
        _points[bufnum[1]].resize(points.size());
        for (int i = 0; i < points.size(); i++) {
            _points[bufnum[1]][i].x() = points[i].x();
            _points[bufnum[1]][i].y() = points[i].y();
            _points[bufnum[1]][i].z() = 0;
        }

        omp_set_lock(&_buflock);
        _bufnum = bufnum;
        _is_reset = true;
        omp_unset_lock(&_buflock);
    }

    void update_buffer(const std::vector<Eigen::Vector2d> &points) {
        std::array<int, 2> bufnum;
        omp_set_lock(&_buflock);
        bufnum = _bufnum;
        omp_unset_lock(&_buflock);

        if ((bufnum[0] != bufnum[1])) {
            return;
        }
        bufnum[1] = (bufnum[1] == 0) ? 1 : 0;
        _points[bufnum[1]].clear();
        _points[bufnum[1]].resize(points.size());
        for (int i = 0; i < points.size(); i++) {
            _points[bufnum[1]][i].x() = points[i].x();
            _points[bufnum[1]][i].y() = points[i].y();
            _points[bufnum[1]][i].z() = 0;
        }

        omp_set_lock(&_buflock);
        _bufnum = bufnum;
        omp_unset_lock(&_buflock);
    }

  protected:
    void set_window() {
        polyscope::options::programName = "point viewer";
        // default 1
        polyscope::options::verbosity = 1;
        // default 1
        polyscope::options::ssaaFactor = 1;
        // default 60
        polyscope::options::maxFPS = 60;
        // default false
        polyscope::options::alwaysRedraw = false;
        // default true
        polyscope::options::buildGui = true;
        // default true
        polyscope::options::openImGuiWindowForUserCallback = true;
        // default true
        polyscope::options::usePrefsFile = true;
    }
    void set_scene() {
        // default false
        polyscope::options::autocenterStructures = false;
        // default false
        polyscope::options::autoscaleStructures = false;
    }
    void set_scene_after() {
        // default true
        polyscope::options::automaticallyComputeSceneExtents = false;
        polyscope::state::lengthScale = 1.;
        polyscope::state::boundingBox =
            std::tuple<glm::vec3, glm::vec3>{{0., 0., 0.}, {1., 1., 1.}};
    }
    void set_camera() {
        // polyscope::view::upDir = polyscope::UpDir::ZUp;
        // Turntable, Free, Planar
        polyscope::view::style = polyscope::NavigateStyle::Planar;
        polyscope::view::bgColor = {179. / 255, 229. / 255, 252. / 255, 1};
    }
    void set_ground_after() {
        // GroundPlaneMode::None, GroundPlaneMode::Tile,
        // GroundPlaneMode::TileReflection, or GroundPlaneMode::ShadowOnly.
        polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
        // defaut 0
        polyscope::options::groundPlaneHeightFactor = 0;
        // defaut .25
        polyscope::options::shadowDarkness = 0.25;
    }
    void set_callback() {
        polyscope::state::userCallback =
            std::function<void()>([&] { callback(); });
    }
    void callback() {
        std::array<int, 2> bufnum;
        bool is_reset;

        omp_set_lock(&_buflock);
        bufnum = _bufnum;
        is_reset = _is_reset;
        omp_unset_lock(&_buflock);

        if (is_reset) {
            if (bufnum[0] != bufnum[1]) {
                polyscope::registerPointCloud(_item_name, _points[bufnum[1]])
                    ->setPointRadius(_radius)
                    ->setPointColor({13. / 255, 71. / 255, 161. / 255});
            }
        } else {
            if (bufnum[0] != bufnum[1]) {
                polyscope::getPointCloud(_item_name)
                    ->updatePointPositions(_points[bufnum[1]]);
            }
        }
        bufnum[0] = bufnum[1];
        is_reset = false;

        omp_set_lock(&_buflock);
        _bufnum = bufnum;
        _is_reset = is_reset;
        omp_unset_lock(&_buflock);
    }

  private:
    std::array<std::vector<Eigen::Vector3d>, 2> _points;
    // displaying buffer, to-display buffer
    std::array<int, 2> _bufnum = {0, 0};
    std::string _item_name = "points";
    bool _is_closed = false;
    bool _is_reset = false;
    double _radius = 0.02;
    int _frame_count = 0;

    omp_lock_t _buflock;
};
