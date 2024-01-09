#include "Eigen/Dense"
#include <iostream>
#include <omp.h>
#include <random>

#include <flip.hpp>
#include <viewer.hpp>

using namespace std;

auto sdfbox(Eigen::Vector2d origin, Eigen::Vector2d range) {
    return [origin, range](Eigen::Vector2d x) {
        x = x - origin;
        x = Eigen::Vector2d(std::abs(x.x()), std::abs(x.y()));
        Eigen::Vector2d q = x - range;
        Eigen::Vector2d q_ =
            Eigen::Vector2d((std::max)(q.x(), 0.0), (std::max)(q.y(), 0.0));
        return q_.norm() + (std::min)((std::max)(q.x(), q.y()), 0.0);
    };
};

auto sdfcircle(Eigen::Vector2d origin, double radius) {
    return [origin, radius](Eigen::Vector2d x) {
        return (x - origin).norm() - radius;
    };
};

int main() {
    omp_set_nested(true);

    Viewer viewer;
    const bool &is_closed = viewer.is_closed();
    FLIPSolver solver(256, 5e-3, 8);
    viewer.set_radius(solver.get_radius());

    solver.add_particle(sdfbox(Eigen::Vector2d(0.2, 0.5), Eigen::Vector2d(0.2, 0.5)));
    solver.set_reset_buffer(
        [&viewer](const std::vector<Eigen::Vector2d> &points) {
            viewer.reset_buffer(points);
        });
    solver.set_update_buffer(
        [&viewer](const std::vector<Eigen::Vector2d> &points) {
            viewer.update_buffer(points);
        });

#pragma omp parallel sections num_threads(2) default(shared)
    {
#pragma omp section
        { viewer.show(); }
#pragma omp section
        { solver.run(is_closed); }
    }
    return 0;
}