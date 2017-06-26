from math import log

def nlogn(n):
    return n*log(n) if n!=0 else 0

def add_count(count, diction, pos):
    if count in diction:
        diction[count].add(pos)
    else:
        diction[count] = set([pos])

class Cluster(object):
    def __init__(self, dim, i, t):
        self.points = []                                #list of points in the cluster
        self.cost = 0.0                                 #the cost of the cluster
        self.dim_counts = [0 for _ in range(dim)]       #the number of 1s in each column
        self.point_count = 0                            #the total number of points
        self.cost_parts = [0.0 for _ in range(dim)]     #the saved parts of the cost that pertain to a certain column
        self.representant = set()                       #the point that represents the cluster
        self.rep_sum = 0                                #the number of 1s in the representant
        self.dis_sum = 0                                #the sum of the differences between the representant and the cluster points
        self.dim = dim                                  #the dimension of the space
        self.number = i
        self.dim_buckets = dict()
        self.t = t

    def __str__(self):
        return str(self.number) + ":" + str(self.point_count) + " " + str(self.dim_counts) + \
            " " + str(self.representant) + " " + str(self.dis_sum) + "\n" + str(self.points) + "\n"

    def addPoint(self, point):                          #adds a point
        self.getCostWith(point, update=True)
        self.points.append(point)
        for pos in point:
            if self.dim_counts[pos] != 0:
                self.dim_buckets[self.dim_counts[pos]].remove(pos)
            self.dim_counts[pos] += 1
            add_count(self.dim_counts[pos], self.dim_buckets, pos)
        self.point_count += 1

    def removePoint(self, point):                       #removes a point
        self.getCostWithout(point, update=True)
        self.points.remove(point)
        for pos in point:
            self.dim_buckets[self.dim_counts[pos]].remove(pos)
            self.dim_counts[pos] -= 1
            add_count(self.dim_counts[pos], self.dim_buckets, pos)
        self.point_count -= 1

    def getCostWith(self, point, update=False):
        """
        Returns the cost that would occur if the point was added
        if update=True, it updates the internal values with the modifications
        observed by the adding of the point
        """
        cost = self.cost
        dis_sum = self.dis_sum
        rep_sum = self.rep_sum
        difference_set = self.representant - point
        for pos in difference_set:
            if self.dim_counts[pos] <= self.t*(self.point_count+1):
                dis_sum += 2*self.dim_counts[pos]-self.point_count
                cost -= nlogn(self.dim_counts[pos]) - self.cost_parts[pos]
                rep_sum -= 1
                if update:
                    self.representant.remove(pos)
                    self.cost_parts[pos] = nlogn(self.dim_counts[pos])
            else:
                dis_sum += 1
                cost -= nlogn(self.point_count - self.dim_counts[pos]+1) - self.cost_parts[pos]
                if update:
                    self.cost_parts[pos] = nlogn(self.point_count - self.dim_counts[pos]+1)
        for pos in point:
            if self.dim_counts[pos]+1 <= self.t*(self.point_count+1):
                dis_sum += 1
                cost -= nlogn(self.dim_counts[pos]+1) - self.cost_parts[pos]
                if update:
                    self.cost_parts[pos] = nlogn(self.dim_counts[pos]+1)
            elif self.dim_counts[pos] <= self.t*self.point_count or self.point_count==0:
                dis_sum += self.point_count-2*self.dim_counts[pos]
                cost -= nlogn(self.point_count-self.dim_counts[pos]) - self.cost_parts[pos]
                rep_sum += 1
                if update:
                    self.representant.add(pos)
                    self.cost_parts[pos] = nlogn(self.point_count-self.dim_counts[pos])
        cost += nlogn(rep_sum) - nlogn(self.rep_sum)
        #cost += (rep_sum - self.rep_sum)*log(self.dim)
        cost += nlogn(dis_sum) - nlogn(self.dis_sum)
        if update:
            self.cost = cost
            self.rep_sum = rep_sum
            self.dis_sum = dis_sum
        return cost

    def getCostWithout(self, point, update=False):
        """
        Returns the cost that would occur if the point was removed
        if update=True, it updates the internal values with the modifications
        observed by the removal of the point
        """
        cost = self.cost
        dis_sum = self.dis_sum
        rep_sum = self.rep_sum
        for pos in self.representant - point:
            dis_sum -= 1
            cost -= nlogn(self.point_count - self.dim_counts[pos] - 1) - self.cost_parts[pos]
            if update:
                self.cost_parts[pos] = nlogn(self.point_count - self.dim_counts[pos] - 1)
        for pos in point:
            if self.dim_counts[pos] <= self.t*self.point_count:
                dis_sum -= 1
                cost -= nlogn(self.dim_counts[pos]-1) - self.cost_parts[pos]
                if update:
                    self.cost_parts[pos] = nlogn(self.dim_counts[pos]-1)
            elif self.dim_counts[pos]-1 <= self.t*(self.point_count-1):
                dis_sum += 2*self.dim_counts[pos]-self.point_count-1
                cost -= nlogn(self.dim_counts[pos]-1) - self.cost_parts[pos]
                rep_sum -= 1
                if update:
                    self.representant.remove(pos)
                    self.cost_parts[pos] = nlogn(self.dim_counts[pos]-1)
        border_val = int(self.point_count*self.t-0.0000001)
        if border_val >= (self.point_count-1)*self.t:
            for i in self.dim_buckets.get(border_val, set()):
                if i not in point:
                    dis_sum += self.point_count-2*self.dim_counts[i]-1
                    rep_sum += 1
                    cost -= nlogn(self.point_count-self.dim_counts[i]-1) - self.cost_parts[i]
                    if update:
                        self.representant.add(i)
                        self.cost_parts[i] = nlogn(self.point_count-self.dim_counts[i]-1)
        cost += nlogn(rep_sum) - nlogn(self.rep_sum)
        #cost += (rep_sum - self.rep_sum)*log(self.dim)
        cost += nlogn(dis_sum) - nlogn(self.dis_sum)
        if update:
            self.cost = cost
            self.dis_sum = dis_sum
            self.rep_sum = rep_sum
        return cost

    def verify(self):
        """
        debug function, checks the integrity of the calculated values
        """
        goods = []
        if not all(goods):
            blad
        goods.append(self.point_count == len(self.points))
        if not all(goods):
            blad
        goods.append(sum(len(point) for point in self.points) == sum(val for val in self.dim_counts))
        if not all(goods):
            blad
        goods.append(len(self.representant) == self.rep_sum)
        if not all(goods):
            blad
        goods.append(all(i in self.representant for i, count in enumerate(self.dim_counts) if count >= self.t*self.point_count))
        if not all(goods):
            blad
        goods.append(all(i not in self.representant for i, count in enumerate(self.dim_counts) if count < self.t*self.point_count))
        if not all(goods):
            blad
        goods.append(self.dis_sum == sum((self.point_count - count) if i in self.representant else count for i, count in enumerate(self.dim_counts)))
        if not all(goods):
            print self, sum((self.point_count - count) if i in self.representant else count for i, count in enumerate(self.dim_counts))
            blad
