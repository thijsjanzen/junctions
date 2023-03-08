//
//  random_functions.hpp
//
//
//  Created by Thijs Janzen on 05/03/2018.
//
//

#ifndef random_functions_hpp
#define random_functions_hpp

#include <random>
#include <vector>
#include <thread>
#include <chrono>
#include <functional>

struct rnd_t {
  std::mt19937 rndgen_;

  rnd_t() {
    std::random_device rd;
    rndgen_ = std::mt19937(rd());
  }

  rnd_t(unsigned int seed) {
    rndgen_ = std::mt19937(seed);
  }

  void set_seed(unsigned int s) {
    rndgen_ = std::mt19937(s);
  }

  std::uniform_real_distribution<> unif_dist =
    std::uniform_real_distribution<>(0, 1.0);

  double uniform() {
    return unif_dist(rndgen_);
  }

  int random_number(unsigned int n) {
    return std::uniform_int_distribution<> (0, n-1)(rndgen_);
  }

  int poisson(double lambda) {
    return std::poisson_distribution<int>(lambda)(rndgen_);
  }

  std::vector<double> generate_random_markers(int number_of_markers) {
    std::vector<double> markers(number_of_markers);
    for(int i = 0; i < number_of_markers; ++i) {
      markers[i] = uniform();
    }
    std::sort(markers.begin(), markers.end());

    return(markers);
  }
};

struct emp_genome {
  std::vector< double > cdf_;
  std::vector< double > pos;

  emp_genome() {
  }

  emp_genome(const emp_genome& other) {
    cdf_ = other.cdf_;
  }

  emp_genome& operator=(const emp_genome& other) {
    if (this != &other) {
      cdf_ = other.cdf_;
    }
    return *this;
  }

  template <typename T>
  emp_genome(const std::vector<T>& positions) {
    pos = positions;
    double total_sum = std::accumulate(positions.begin(),
                                       positions.end(), 0.0);
    double s = 0.0;
    double mult = 1.0 / total_sum;
    cdf_.resize(positions.size());
    for (size_t i = 0; i < positions.size(); ++i) {
      s += positions[i] * mult;
      cdf_[i] = s;
    }
    return;
  }

  size_t index_from_cdf(double p) const {
    // find index belonging to p
    return static_cast<size_t>(std::distance(cdf_.begin(),
                                             std::lower_bound(cdf_.begin(),
                                                              cdf_.end(),
                                                              p)));
  }

  std::vector< size_t > recompos(double morgan,
                                 rnd_t& rndgen) const {
    size_t num_break_points = rndgen.poisson(morgan);
    std::vector< size_t > indices;
    for(size_t i = 0; i < num_break_points; ++i) {
      auto found_index = index_from_cdf(rndgen.uniform());
      if (found_index > 0) {
        indices.push_back(found_index);
      }
    }
    std::sort(indices.begin(), indices.end());
    indices.push_back(cdf_.size());
    return indices;
  }
};

inline int get_seed() {
  const auto tt = static_cast<int64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  auto tid = std::this_thread::get_id();
  const uint64_t e3{ std::hash<std::remove_const_t<decltype(tid)>>()(tid) };
  return static_cast<int>(tt + e3);
}


#endif /* random_functions_hpp */
