

#include "cluster.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <list>
#include <iostream>

using namespace std;

double Cluster::getCost() {
    return cost;
}

Cluster::Cluster(int d, double t, int n):
    points(list< Point*>()), cost(0), dim_counts(new int[d]),
    point_count(0), cost_parts(new double[d]), representative(Point()),
    rep_sum(0), dis_sum(0), dim(d), dim_buckets(new set<int>[n]), t(t){
    for (int i=0;i<d;++i) {
        cost_parts[i] = 0.0;
        dim_counts[i] = 0;
    }
}

//utility function
static double nlogn(int n) {
    return n==0 ? 0 : log(n)*n;
}

void Cluster::addPoint( Point* point) {
    getCostWith(point, true);
    points.push_back(point);
    for (Point::iterator it = point->begin(); it != point->end(); ++it) {
        if (dim_counts[*it] != 0) {
            dim_buckets[dim_counts[*it]].erase(*it);
        }
        ++dim_counts[*it];
        dim_buckets[dim_counts[*it]].insert(*it);
    }
    ++point_count;
}

void Cluster::removePoint( Point* point) {
    getCostWithout(point, true);
    points.remove(point);
    for (Point::iterator it = point->begin(); it != point->end(); ++it) {
        dim_buckets[dim_counts[*it]].erase(*it);
        --dim_counts[*it];
        if (dim_counts[*it] != 0) {
            dim_buckets[dim_counts[*it]].insert(*it);
        }
    }
    --point_count;
}

//utility function, removes the first occurence of a value in a vector
static void removeFromVector(Point& point, int val) {
    for (Point::iterator it = point.begin(); it != point.end(); ++it) {
        if (*it == val) {
            point.erase(it);
            break;
        }
    }
}

double Cluster::getCostWith( Point* point, bool update) {
    double cost_ = this->cost;              //local copies of the cluster variables
    int dis_sum_ = this->dis_sum;
    int rep_sum_ = this->rep_sum;
    vector<int> difference_set(representative.size());

    vector<int>::iterator itend = set_difference(representative.begin(),
        representative.end(), point->begin(), point->end(), difference_set.begin());
    for (vector<int>::iterator it = difference_set.begin(); it != itend; ++it) {
        //position is in the representative but not in the point being added
        if (dim_counts[*it] < t*(point_count+1)) {
            //new point causes the position to be removed from the representative
            dis_sum_ += 2*dim_counts[*it]-point_count;
            cost_ -= nlogn(dim_counts[*it]) - cost_parts[*it];
            --rep_sum_;
            if (update) {
                removeFromVector(representative, *it);
                cost_parts[*it] = nlogn(dim_counts[*it]);
            }
        } else {
            //new point does not cause a change in the representative
            ++dis_sum_;
            cost_ -= nlogn(point_count - dim_counts[*it]+1) - cost_parts[*it];
            if (update)
                cost_parts[*it] = nlogn(point_count - dim_counts[*it]+1);
        }
    }

    for (Point::iterator it = point->begin(); it != point->end(); ++it) {
        if (dim_counts[*it]+1 < t*(point_count+1)) {
            //the position is not in the representative but is in the new point
            ++dis_sum_;
            cost_ -= nlogn(dim_counts[*it]+1) - cost_parts[*it];
            if (update)
                cost_parts[*it] = nlogn(dim_counts[*it]+1);
        } else if (dim_counts[*it] < t*point_count || point_count==0) {
            //the position is in the new point, which causes it
            //to cross the treshold and be added to the representative
            dis_sum_ += point_count-2*dim_counts[*it];
            cost_ -= nlogn(point_count-dim_counts[*it]) - cost_parts[*it];
            ++rep_sum_;
            if (update) {
                representative.push_back(*it);
                cost_parts[*it] = nlogn(point_count-dim_counts[*it]);
            }
        }
    }

    cost_ += nlogn(rep_sum_) - nlogn(rep_sum);
    cost_ += nlogn(dis_sum_) - nlogn(dis_sum);
    if (update) {
        cost = cost_;
        rep_sum = rep_sum_;
        dis_sum = dis_sum_;
    }
    return cost_;
}

//utility function, returns whether a value is present in a vector
static bool valInPoint( Point* point, int val) {
    for (Point::iterator it = point->begin(); it != point->end(); ++it) {
        if (*it == val) {
            return true;
        }
    }
    return false;
}

double Cluster::getCostWithout( Point* point, bool update) {
    double cost_ = this->cost;              //local copies of the cluster variables
    int dis_sum_ = this->dis_sum;
    int rep_sum_ = this->rep_sum;
    vector<int> difference_set(representative.size());

    vector<int>::iterator itend = set_difference(representative.begin(),
        representative.end(), point->begin(), point->end(), difference_set.begin());
    for (vector<int>::iterator it = difference_set.begin(); it != itend; ++it) {
        //position is in the representative but not in the point being removed
        --dis_sum_;
        cost_ -= nlogn(point_count - dim_counts[*it] - 1) - cost_parts[*it];
        if (update)
               cost_parts[*it] = nlogn(point_count - dim_counts[*it] - 1);
    }

    for (Point::iterator it = point->begin(); it != point->end(); ++it) {
        if (dim_counts[*it] < t*point_count) {
            //position is in the removed point but does not change the representative
            --dis_sum_;
            cost_ -= nlogn(dim_counts[*it]-1) - cost_parts[*it];
            if (update)
                cost_parts[*it] = nlogn(dim_counts[*it]-1);
        } else if (dim_counts[*it]-1 < t*(point_count-1)) {
            //position in the removed point causes it to be removed from the representative
            dis_sum_ += 2*dim_counts[*it] - point_count - 1;
            cost_ -= nlogn(dim_counts[*it] - 1) - cost_parts[*it];
            --rep_sum_;
            if (update) {
                removeFromVector(representative, *it);
                cost_parts[*it] = nlogn(dim_counts[*it] - 1);
            }
        }
    }

    //special case: checking if removing the point would
    //cause a position to be added to the representative
    int border_val = (int)(point_count*t-0.0000001);
    if (border_val >= (point_count-1)*t) {
        set<int>& bucket = dim_buckets[border_val];
        for (set<int>::iterator it = bucket.begin(); it != bucket.end(); ++it)
            if (valInPoint(point, *it)) {
                dis_sum_ += point_count-2*dim_counts[*it]-1;
                ++rep_sum_;
                cost_ -= nlogn(point_count-dim_counts[*it]-1) - cost_parts[*it];
                if (update) {
                    representative.push_back(*it);
                    cost_parts[*it] = nlogn(point_count-dim_counts[*it]-1);
                }
            }
    }

    cost_ += nlogn(rep_sum_) - nlogn(rep_sum);
    cost_ += nlogn(dis_sum_) - nlogn(dis_sum);
    if (update) {
        cost = cost_;
        rep_sum = rep_sum_;
        dis_sum = dis_sum_;
    }
    return cost_;
}

