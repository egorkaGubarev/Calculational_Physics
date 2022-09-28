#include <iostream>
#include <cmath>
#include <fstream>
#include <list>
#include <vector>

typedef unsigned int uint;
typedef float my_type;

template <typename type>
struct Point
{
    type x_, y_, z_;

    Point(const type x, const type y, const type z) : x_{ x }, y_{ y }, z_{ z } {}

    std::ofstream& Log(std::ofstream& out) const
    {
        out << x_ << ' ' << y_ << ' ' << z_ << '\n';
        return out;
    }

    ~Point() = default;
};

template <typename type>
struct Segment
{
    const Point<type> begin_, end_;

    Segment(const Point<type> begin, const Point<type> end): begin_{begin}, end_{end} {}

    ~Segment() = default;
};

template <typename type>
struct Node
{
    bool is_rigid_;
    Point<type> point_;
    std::vector<std::reference_wrapper<Node<type>>> connected_nodes_;

    Node(Point<type> point, const bool is_rigid) : point_{ point }, is_rigid_{ is_rigid } {}

    void connect(Node<type>& other)
    {
        connected_nodes_.push_back(other);
        std::vector<std::reference_wrapper<Node<type>>>& connections = other.connected_nodes_;
        connections.push_back(*this);
    }
};

template <typename type>
struct Triangle
{
    const Segment<type> segm_1_, segm_2_, segm_3_;

    Triangle(const Segment<type> segm_1, const Segment<type> segm_2, const Segment<type> segm_3):
    segm_1_{segm_1}, segm_2_{segm_2}, segm_3_{segm_3} {}

    ~Triangle() = default;
};

template <typename type>
struct Circuit
{
    my_type x_max_, x_min_, y_max_, y_min_, z_max_, z_min_, x_center_, y_center_, z_center_;
    std::list<Node<type>> soft_nodes_, rigid_nodes_;

    Circuit(Point<type> generator(const type t), const uint dots)
    {
        const type step = (type) 1 / dots;
        Point<type> zero = generator(0);
        create_rigid_node(zero);

        for (type t = step; t < 1; t += step) {
            Point<type> new_point = generator(t);
            create_rigid_node(new_point);
            Node<type>& prev = get_last_but_one_rigid_node();
            Node<type>& current = get_last_rigid_node();
            current.connect(prev);
        }

        connect_last_and_zero();
    }

    void get_borders_coords()
    {
        bool x_max_found = false, x_min_found = false, y_max_found = false, y_min_found = false; bool z_min_found = false, z_max_found = false;

        for (const Node<my_type>& node : rigid_nodes_) {
            const Point<my_type>& point = node.point_;

            const my_type x = point.x_;
            const my_type y = point.y_;
            const my_type z = point.z_;

            if (!x_max_found || x > x_max_) {
                x_max_ = x;
                x_max_found = true;
            }

            if (!x_min_found || x < x_min_) {
                x_min_ = x;
                x_min_found = true;
            }

            if (!y_max_found || y > y_max_) {
                y_max_ = y;
                y_max_found = true;
            }

            if (!y_min_found || y < y_min_) {
                y_min_ = y;
                y_min_found = true;
            }

            if (!z_max_found || z > z_max_) {
                z_max_ = z;
                z_max_found = true;
            }

            if (!z_min_found || z < z_min_) {
                z_min_ = z;
                z_min_found = true;
            }
        }
    }

    void get_center()
    {
        x_center_ = (x_min_ + x_max_) / 2;
        y_center_ = (y_min_ + y_max_) / 2;
        z_center_ = (z_min_ + z_max_) / 2;
    }

    void create_soft_node(const Point<type> point)
    {
        Node<type> node(point, false);
        soft_nodes_.push_back(node);
    }

    void create_rigid_node(const Point<type> point)
    {
        Node<type> node(point, true);
        rigid_nodes_.push_back(node);
    }

    void create_primary_node()
    {
        Point<type> center(x_center_, y_center_, z_center_);
        create_soft_node(center);
    }

    // Use when only primary node is created
    void create_primary_surf()
    {
        typename std::list<Node<type>>::iterator begin = std::begin(soft_nodes_);
        Node<type>& center = *begin;
        std::vector<std::reference_wrapper<Node<type>>>& connected_nodes = center.connected_nodes_;

        for (Node<type>& to_connect : rigid_nodes_) {
            connected_nodes.push_back(to_connect);
            std::vector<std::reference_wrapper<Node<type>>>& reverse_connect = to_connect.connected_nodes_;
            reverse_connect.push_back(center);
        }
    }

    // Not ready
    void collect_triangles()
    {
        for (const Node<type>& point: soft_nodes_) {
            const std::vector<std::reference_wrapper<Node<type>>>& second_points = point.connected_nodes_;

        }
    }

    void get_surf()
    {
        get_borders_coords();
        get_center();
        create_primary_node();
        create_primary_surf();
    }

    Node<type>& get_last_rigid_node()
    {
        typename std::list<Node<type>>::iterator end = std::end(rigid_nodes_);
        typename std::list<Node<type>>::iterator last = std::prev(end);
        Node<type>& node = *last;
        return node;
    }

    // Use when you have at least two nodes only
    Node<type>& get_last_but_one_rigid_node()
    {
        typename std::list<Node<type>>::iterator end = std::end(rigid_nodes_);
        typename std::list<Node<type>>::iterator last = std::prev(end, 2);
        Node<type>& node = *last;
        return node;
    }

    void connect_last_and_zero()
    {
        Node<type>& last = get_last_rigid_node();
        typename std::list<Node<type>>::iterator begin = std::begin(rigid_nodes_);
        Node<type>& first = *begin;
        last.connect(first);
    }

    std::ofstream& log_rigid(std::ofstream& out) const
    {
        for (const Node<type>& node : rigid_nodes_) {
            const Point<type>& point = node.point_;

            const type x = point.x_;
            const type y = point.y_;
            const type z = point.z_;

            out << x << ' ' << y << ' ' << z << '\n';
        }

        return out;
    }

    std::ofstream& log_soft(std::ofstream& out) const
    {
        for (const Node<type>& node : soft_nodes_) {
            const Point<type>& point = node.point_;

            const type x = point.x_;
            const type y = point.y_;
            const type z = point.z_;

            out << x << ' ' << y << ' ' << z << '\n';
        }

        return out;
    }

    std::ofstream& log_connect(std::ofstream& out) const
    {
        for (const Node<type>& node : soft_nodes_) {
            const Point<type>& point = node.point_;
            const std::vector<std::reference_wrapper<Node<type>>>& connect = node.connected_nodes_;

            const type x = point.x_;
            const type y = point.y_;
            const type z = point.z_;

            for (const Node<type>& connect_node : connect) {
                const Point<type>& connect_point = connect_node.point_;

                const type connect_x = connect_point.x_;
                const type connect_y = connect_point.y_;
                const type connect_z = connect_point.z_;

                out << x << ' ' << y << ' ' << z << '\n';
                out << connect_x << ' ' << connect_y << ' ' << connect_z << '\n';
            }
        }

        return out;
    }

    std::ofstream& log(std::ofstream& out) const
    {
        log_rigid(out);
        log_soft(out);
        log_connect(out);

        return out;
    }

    ~Circuit() = default;
};

template <typename type>
Point<type> circle(const type t)
{
    const type pi = (type) std::acos(-1);
    const type angle = 2 * pi * t;

    const type x = std::cos(angle);
    const type y = std::sin(angle);
    const type z = 0;

    const Point<type> point(x, y, z);
    return point;
}

template <typename type>
Point<type> deform_circle(const type t)
{
    const type vertic_amplit = (type) 0.1;
    const type vertic_frequency = 10;

    const type pi = (type)std::acos(-1);
    const type angle = 2 * pi * t;

    const type x = std::cos(angle);
    const type y = std::sin(angle);
    const type z = vertic_amplit * std::sin(angle * vertic_frequency);

    const Point<type> point(x, y, z);
    return point;
}

template <typename type>
Point<type> eight(const type t)
{
    const type vertic_amplit = (type)0.1;
    const type vertic_frequency = 0.25;

    const type pi = (type)std::acos(-1);
    const type angle = 4 * pi * t;
    type y = 0;
    const type x = (type)std::sin(angle);

    if (angle < 2 * pi) {
        y = (type) (1 - std::cos(angle));
    }
    else {
        y = (type) (std::cos(angle) - 1);
    }

    const type z = vertic_amplit * std::sin(angle * vertic_frequency);
    const Point<type> point(x, y, z);
    return point;
}

int main()
{
    const uint dots_on_circ = 100;
    std::ofstream aux, file_circle, file_deform_circle, file_eight;

    aux.open("aux_info.txt");
    file_circle.open("circle.txt");
    file_deform_circle.open("deform_circle.txt");
    file_eight.open("eight.txt");

    Circuit<my_type> circ_circle(circle, dots_on_circ), circ_deform_circle(deform_circle, dots_on_circ), circ_eight(eight, dots_on_circ);

    aux << dots_on_circ << '\n';
    aux.close();
    std::vector<std::reference_wrapper<std::ofstream>> files{file_circle, file_deform_circle, file_eight};
    std::vector<std::reference_wrapper<Circuit<my_type>>> circs{ circ_circle, circ_deform_circle, circ_eight };
    const size_t circs_amount = std::size(circs);

    for (size_t circ_numb = 0; circ_numb < circs_amount; ++circ_numb) {
        Circuit<my_type>& circ = circs[circ_numb];
        circ.get_surf();
        std::ofstream& file = files[circ_numb];
        circ.log(file);
    }
    
    for (std::ofstream& file : files) {
        file.close();
    }

    return 0;
}