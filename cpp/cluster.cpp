

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
	point_count(0), cost_parts(new double[d]), representant(Point()),
	rep_sum(0), dis_sum(0), dim(d), dim_buckets(new set<int>[n]), t(t){
	for (int i=0;i<d;++i) {
		cost_parts[i] = 0.0;
		dim_counts[i] = 0;
	}
}

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

static void removeFromVector(Point& point, int val) {
	for (Point::iterator it = point.begin(); it != point.end(); ++it) {
		if (*it == val) {
			point.erase(it);
			break;
		}
	}
}

double Cluster::getCostWith( Point* point, bool update) {
	double cost_ = this->cost;
	int dis_sum_ = this->dis_sum;
	int rep_sum_ = this->rep_sum;
	vector<int> difference_set(representant.size());

	vector<int>::iterator itend = set_difference(representant.begin(), representant.end(), point->begin(), point->end(), difference_set.begin());
	for (vector<int>::iterator it = difference_set.begin(); it != itend; ++it) {
		if (dim_counts[*it] < t*(point_count+1)) {
        	dis_sum_ += 2*dim_counts[*it]-point_count;
            cost_ -= nlogn(dim_counts[*it]) - cost_parts[*it];
        	--rep_sum_;
        	if (update) {
				removeFromVector(representant, *it);
            	cost_parts[*it] = nlogn(dim_counts[*it]);
			}
		} else {
			++dis_sum_;
			cost_ -= nlogn(point_count - dim_counts[*it]+1) - cost_parts[*it];
			if (update)
                cost_parts[*it] = nlogn(point_count - dim_counts[*it]+1);
		}
	}
	
	for (Point::iterator it = point->begin(); it != point->end(); ++it) {
		if (dim_counts[*it]+1 < t*(point_count+1)) {
            ++dis_sum_;
            cost_ -= nlogn(dim_counts[*it]+1) - cost_parts[*it];
            if (update)
                cost_parts[*it] = nlogn(dim_counts[*it]+1);
        } else if (dim_counts[*it] < t*point_count || point_count==0) {
            dis_sum_ += point_count-2*dim_counts[*it];
            cost_ -= nlogn(point_count-dim_counts[*it]) - cost_parts[*it];
            ++rep_sum_;
            if (update) {
                representant.push_back(*it);
                cost_parts[*it] = nlogn(point_count-dim_counts[*it]);
			}
		}
	}
	
	cost_ += nlogn(rep_sum_) - nlogn(rep_sum);
	double prevcost = cost_;
	cost_ += nlogn(dis_sum_) - nlogn(dis_sum);
	if (update) {
		cost = cost_;
		rep_sum = rep_sum_;
		dis_sum = dis_sum_;
	}
	return cost_;
}

static bool valInPoint( Point* point, int val) {
	for (Point::iterator it = point->begin(); it != point->end(); ++it) {
		if (*it == val) {
			return true;
		}
	}
	return false;
}

double Cluster::getCostWithout( Point* point, bool update) {
	double cost_ = this->cost;
	int dis_sum_ = this->dis_sum;
	int rep_sum_ = this->rep_sum;
	vector<int> difference_set(representant.size());

	vector<int>::iterator itend = set_difference(representant.begin(), representant.end(), point->begin(), point->end(), difference_set.begin());
	for (vector<int>::iterator it = difference_set.begin(); it != itend; ++it) {
    	--dis_sum_;
        cost_ -= nlogn(point_count - dim_counts[*it] - 1) - cost_parts[*it];
        if (update)
           	cost_parts[*it] = nlogn(point_count - dim_counts[*it] - 1);
	}

	for (Point::iterator it = point->begin(); it != point->end(); ++it) {
		if (dim_counts[*it] < t*point_count) {
            --dis_sum_;
            cost_ -= nlogn(dim_counts[*it]-1) - cost_parts[*it];
            if (update)
                cost_parts[*it] = nlogn(dim_counts[*it]-1);
        } else if (dim_counts[*it]-1 < t*(point_count-1)) {
            dis_sum_ += 2*dim_counts[*it] - point_count - 1;
            cost_ -= nlogn(dim_counts[*it] - 1) - cost_parts[*it];
            --rep_sum_;
			if (update) {
				removeFromVector(representant, *it);
                cost_parts[*it] = nlogn(dim_counts[*it] - 1);
			}
		}
	}
	
	int border_val = (int)(point_count*t-0.0000001);
	if (border_val >= (point_count-1)*t) {
		set<int>& bucket = dim_buckets[border_val];
        for (set<int>::iterator it = bucket.begin(); it != bucket.end(); ++it)
			if (valInPoint(point, *it)) {
                dis_sum_ += point_count-2*dim_counts[*it]-1;
                ++rep_sum_;
                cost_ -= nlogn(point_count-dim_counts[*it]-1) - cost_parts[*it];
                if (update) {
                    representant.push_back(*it);
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

