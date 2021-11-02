/**
 * JStat is a collection of statistical methods
 * */
var jstat = {
	quicksort: function(data) {
		var smaller = [];
		var larger = [];
		var pivot = Math.floor(data.length / 2);
		if (data.length <= 1) return data;
		for (var i = 0; i < data.length; i++) {
			if (i == pivot) {
				continue
			}
			if (data[pivot] > data[i]) {
				smaller.push(data[i]);
			}
			if (data[pivot] < data[i]) {
				larger.push(data[i]);
			}
		}

		var left = this.quicksort(smaller);
		var right = this.quicksort(larger);
		return left.concat(data[pivot]).concat(right);
	},

	/**
	 * Sample min
	 */
	sample_min: function(data) {
		var min = undefined;
		for (var i = 0; i < data.length; i++) {
			if (min == undefined) {
				min = data[i];
			} else {
				if (min > data[i]) {
					min = data[i];
				}
			}
		}
		return min;
	},

	/**
	 * Sample max
	 */
	sample_max: function(data) {
		var max = undefined;
		for (var i = 0; i < data.length; i++) {
			if (max == undefined) {
				max = data[i];
			} else {
				if (max < data[i]) {
					max = data[i];
				}
			}
		}
		return max;
	},

	/**
	 * Computes sample mean
	 */
	sample_mean: function(data) {
		var sum = 0;
		var n = data.length;
		for (var i = 0; i < n; i++) {
			sum += data[i];
		}
		return sum / n;
	},

	/** 
	 * 
	 * Computes sample variance
	 */
	sample_var: function(data) {
		if (data.length <= 1) {
			return 0;
		}

		var mean = this.sample_mean(data);
		var sum = 0;
		var n = data.length - 1;
		for (var i = 0; i < n + 1; i++) {
			sum += (data[i] - mean) * (data[i] - 1);
		}

		return sum / n;
	},

	/**
	 * Computes sample standard deviation
	 */
	sample_std: function(data) {
		return Math.sqrt(this.sample_var(data));
	},

	/**
	 * Compute median
	 */
	median: function(data) {
		if (data.length == 0) {
			return undefined;
		}
		var sorted = this.quicksort(data);
		console.log(sorted);
		if (sorted.length % 2 == 1) {
			return sorted[Math.floor(sorted.length / 2)];
		} else {
			return (sorted[sorted.length / 2 - 1] + sorted[sorted.length / 2]) / 2;
		}
	},

	/**
	 * Returns random value from the uniform distribution
	 */
	runif: function(a, b) {
		return Math.random()*(b - a) + a;
	},

	/**
	 * Quantile function of the uniform distribution
	 */
	quinf: function(a, b, p) {
		return (b - a) * p + a;
	},

	/**
	 * Cumulative distribution function
	 */
	punif: function(a, b, x) {
		if (x < a) return 0;
		if (x > b) return 1;
		return (x - a) / (b - a);
	},

	/**
	 * Kruskal-Wallis test
	 * Time complexity of the test is O(nlogn)
	*/
	kruska_F_score: function(data) {
		/*
		Space complexity: O(nlogn)
		Time complexity: O(nlogn)
		*/
		function quicksort(arr, groups) {
			if (arr == undefined) return [[], []];
			if (arr.length == 1) return [arr, groups];
			if (arr.length == 0) return [arr, groups];
			var midpoint = Math.floor(arr.length / 2);
			var smaller = [];
			var larger = [];
			var smaller_g = [];
			var larger_g = [];
			for (var i = 0; i < arr.length; i++) {
				if (i == midpoint) continue;
				if (arr[midpoint] >= arr[i]) {
					smaller.push(arr[i]);
					smaller_g.push(groups[i])
				}
				if (arr[midpoint] < arr[i]) {
					larger.push(arr[i]);
					larger_g.push(groups[i])
				}
			}
			var [left, left_g] = quicksort(smaller, smaller_g);
			var [right, right_g] = quicksort(larger, larger_g);
			return [left.concat([arr[midpoint]]).concat(right),
				left_g.concat([groups[midpoint]]).concat(right_g)];
		}

		/*
		Space complexity: O(n)
		Time complexity (worst case): O(n)
		*/
		function rank(x) {
			var r = [];
			var avg_rank = 0;
			var avg_rank_count = 0;
			for (var i = 0; i < x.length; i++) {
				avg_rank += (i + 1);
				avg_rank_count += 1;
				if (i + 1 < x.length) {
					if (x[i] != x[i + 1]) {
						for (var j = 0; j < avg_rank_count; j++) {
							r.push(avg_rank / avg_rank_count);
						}
						avg_rank = 0;
						avg_rank_count = 0;
					}
				} else {
					for (var j = 0; j < avg_rank_count; j++) {
						r.push(avg_rank / avg_rank_count);
					}
				}
			}
			return r;
		}

		/*
		Space complexity: O(n)
		Time complexity (worst case): O(nm+nlogn+n+n+m) = O(max(nlogn, nm))
		*/
		function kruskal(arrs) {
			var x = [], r = [], g = [], cr = [];
			//O(nm)
			for (var i = 0; i < arrs.length; i++) {
				x = x.concat(arrs[i]);
				cr.push(0);
				for (var j = 0; j < arrs[i].length; j++) {
					g.push(i+1);
				}
			}
			//O(nlogn)
			var [x, g] = quicksort(x, g);
			//O(n)
			r = rank(x);
			//O(n)
			for (var i = 0; i < x.length; i++) {
				cr[g[i] - 1] += r[i];
			}
			var sum = 0, n = x.length;
			//O(m)
			for (var i = 0; i < arrs.length; i++) {
				sum += cr[i]*cr[i] / arrs[i].length;
			}
			return 12/(n*(n+1))*sum - 3*(n+1);
		}

		return kruskal(data);
	}
}
