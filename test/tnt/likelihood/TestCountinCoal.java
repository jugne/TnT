package tnt.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math.util.MathUtils;
import org.junit.Test;

import beast.util.Randomizer;

public class TestCountinCoal {
	static double tau = 2.0;
	static double Ne = 1.0;
	static double iterations = 100000000;

	@Test
	public void CountSubCoal() {

		// tips are labeled 1..5, direction doesn't matter: (1,2)=(2,1)
		// # of coalescences, where 5 lineages coalesce to 2 (5->2) in any order
		int total_ways_to_6_3 = 0;

		// # of coalescences, where (6->3) and there are three coalescent "clusters"
		// formed: (1,2,3) and (4,5) and (6). Coalescent order within the clusters is
		// not
		// important.
		int ways_to_123_and_45_and_6_in6_any = 0;

		// # of coalescences, where (6->2) and there are two coalescent "clusters"
		// formed: (1,2,3) and (4,5) and (6). Coalescent order within the clusters is
		// case 0: ((1,2),3) and (4,5) and (6)
		// case 1: ((1,3),2) and (4,5) and (6)
		// case 2: ((2,3),1) and (4,5) and (6)
		// Coalescent order between the clusters is ANY.
		int ways_to_123_and_45_and_6_in6_ordered_0 = 0;
		int ways_to_123_and_46_and_6_in6_ordered_1 = 0;
		int ways_to_123_and_46_and_6_in6_ordered_2 = 0;

		// As above case 0, but order between the clusters is SPECIFIC
		int ways_to_123_and_45_and_6_in6_ordered_00 = 0;

		// tips are labeled 1..6, direction doesn't matter: (1,2)=(2,1)
		// # of coalescences, where 6 lineages coalesce to 2 (6->2) in any order
		int total_ways_to_6_2 = 0;

		// # of coalescences, where (6->2) and there are two coalescent "clusters"
		// formed: (1,2,3) and (4,5, 6). Coalescent order within the clusters is not
		// important.
		int ways_to_123_and_456_in6_any = 0;

		for (int j = 0; j < iterations; j++) {
			int n = 6;
		int[] ind = new int[n];
			for (int i = 1; i <= n; i++) {
			ind[i-1]=i;
		}
		

		List<Integer[]> coal = new ArrayList<Integer[]>();
			List<Integer[]> coal_list = new ArrayList<Integer[]>();
		
			for (int i = 1; i <= n; i++) {
			coal.add(new Integer[] { i });
		}
		
		double duplicateTime = 0;
			double stopTime = 0 + tau;

		while (duplicateTime < stopTime) {
			int nrLineages = n;
			double deltaT = Randomizer
					.nextExponential(nrLineages * (nrLineages - 1) * 0.5 * 1 / Ne);
			if (duplicateTime + deltaT < stopTime) {
				int k = nrLineages;
				Integer[] node1 = coal.get(Randomizer.nextInt(k));
				Integer[] node2;
				do {
					node2 = coal.get(Randomizer.nextInt(k));
				} while (node2 == node1);

					Integer[] c = concat(node1, node2);
					coal_list.add(concat(node1, node2));
					coal.add(c);
					coal.remove(node1);
					coal.remove(node2);
				
				n = coal.size();
		}
			duplicateTime += deltaT;
		}
		
			List<Integer> lengths = new ArrayList<>();
			if (coal_list.size() == 3) {
				for (int i = 0; i < coal_list.size(); i++) {
					lengths.add(coal_list.get(i).length);
				}
				if (lengths.contains(3) && lengths.contains(2) && !lengths.contains(4)) {
					total_ways_to_6_3 += 1;
				for (int i = 0; i < coal_list.size(); i++) {
						List<Integer> list = Arrays.asList(coal_list.get(i));
						if (list.size() == 3 && list.contains(1) && list.contains(2) && list.contains(3)) {
							for (int k = 0; k < coal_list.size(); k++) {
								if (k != i) {
								List<Integer> list2 = Arrays.asList(coal_list.get(k));
								if (list2.size() == 2 && list2.contains(5) && list2.contains(4)) {
										ways_to_123_and_45_and_6_in6_any += 1;
										for (int l = 0; l < coal_list.size(); l++) {
											if (l != i && l != k) {
												List<Integer> list3 = Arrays.asList(coal_list.get(l));
												if (list3.size() == 2 && list3.contains(1) && list3.contains(2)) {
													ways_to_123_and_45_and_6_in6_ordered_0 += 1;
													if (l < i && l < k && k < i)
														ways_to_123_and_45_and_6_in6_ordered_00 += 1;
//												System.out.println(i + "," + k + "," + l);
												}
												if (list3.size() == 2 && list3.contains(1) && list3.contains(3)) {
													ways_to_123_and_46_and_6_in6_ordered_1 += 1;
//												System.out.println(i + "," + k + "," + l);
												}
												if (list3.size() == 2 && list3.contains(2) && list3.contains(3)) {
													ways_to_123_and_46_and_6_in6_ordered_2 += 1;
//												System.out.println(i + "," + k + "," + l);
												}
											}
										}
									}
								}
							}
						}
					}
				}
				if (lengths.contains(4))
					total_ways_to_6_3 += 1;
				if (lengths.contains(2) && !lengths.contains(3) && !lengths.contains(1) && !lengths.contains(4))
					total_ways_to_6_3 += 1;
			} else if (coal_list.size() == 4) {
				for (int i = 0; i < coal_list.size(); i++) {
					lengths.add(coal_list.get(i).length);
				}
				if (lengths.contains(3) && !lengths.contains(4) && !lengths.contains(5)) {
					total_ways_to_6_2 += 1;
					for (int i = 0; i < coal_list.size(); i++) {
						List<Integer> list = Arrays.asList(coal_list.get(i));
						if (list.size() == 3 && list.contains(1) && list.contains(2) && list.contains(3)) {
							for (int k = 0; k < coal_list.size(); k++) {
								if (k != i) {
									List<Integer> list2 = Arrays.asList(coal_list.get(k));
									if (list2.size() == 3 && list2.contains(4) && list2.contains(5)
											&& list2.contains(6)) {
										ways_to_123_and_456_in6_any += 1;
									}
								}
							}
						}
					}
				}
				if (lengths.contains(4) && !lengths.contains(5))
					total_ways_to_6_2 += 1;
				if (lengths.contains(5))
					total_ways_to_6_2 += 1;

			}

		}
		System.out.println("Total ways 6->3: ");
		System.out.println("Empirical: " + total_ways_to_6_3 / iterations +
				"\nAnalytical: " + gUp(6, 3, tau, Ne));

		System.out.println();
		System.out.println("Total ways 6->3, ANY order within sub-coalescent clusters ((1,2),3) and (4,5): ");
		System.out.println("Empirical: " + ways_to_123_and_45_and_6_in6_any / iterations +
				"\nAnalytical: " + coefficient(Arrays.asList(1, 2, 3), 6, 3));

		System.out.println();
		System.out.println("Total ways 6->3, SPECIFIC order within sub-coalescent clusters ((1,2),3) and (4,5),"
				+ "\n ANY order of coalescence between sub-clusters: ");
		System.out.println("Empirical: " + ways_to_123_and_45_and_6_in6_ordered_0 / iterations +
				"\nAnalytical: " + gUp(6, 3, tau, Ne) * (1.0 / waysToCoal(6, 3)) * binomialInt(3 - 1 + 2 - 1, 3 - 1));

		System.out.println();
		System.out.println("Total ways 6->3, specific order within sub-coalescent clusters ((1,2),3) and (4,5),"
				+ "\n SPECIFIC order of coalescence between sub-clusters: ");
		System.out.println("Empirical: " + ways_to_123_and_45_and_6_in6_ordered_00 / iterations +
				"\nAnalytical: " + gUp(6, 3, tau, Ne) * (1.0 / waysToCoal(6, 3)));

		System.out.println();
		System.out.println("---------------------------------------------------------------------------");
		System.out.println("Total ways 6->2: ");
		System.out.println("Empirical: " + total_ways_to_6_2 / iterations +
				"\nAnalytical: " + gUp(6, 2, tau, Ne));
		System.out.println("Total ways 6->2, ANY order within sub-coalescent clusters ((1,2),3) and ((4,5),6): ");
		System.out.println("Empirical: " + ways_to_123_and_456_in6_any / iterations +
				"\nAnalytical: " + coefficient(Arrays.asList(1, 3, 3), 6, 2));

	}

	public static <T> T[] concat(T[] first, T[] second) {
		T[] result = Arrays.copyOf(first, first.length + second.length);
		System.arraycopy(second, 0, result, first.length, second.length);
		return result;
	}

	// Calculates (g_{i,j}(tau)/I_{i,j})*
	// I_{k_1,1}*I_{k_2,1}*binom(k_1-1+k_2-1,k_1-1)
	static double coefficient(List<Integer> k, int l_1, int l_2) {
		double ans = 0.0;
		double mult = 1.0;
		int sum = k.stream().mapToInt(Integer::intValue).sum();
		for (int s = 0; s < k.size(); s++) {
			mult *= waysToCoal(k.get(s), 1);
			// W factor from NOAH A. ROSENBERG
			mult *= binomialInt(sum - k.size(), k.get(s) - 1);
			sum -= (k.get(s) - 1);
		}
		ans += (1.0 / waysToCoal(l_1, l_2)) * mult
				* gUp(l_1, l_2, tau, Ne);

		return ans;
	}

	static double waysToCoal(int i, int j) {
		double ans = MathUtils.factorial(i);
		ans *= MathUtils.factorial(i - 1);
		ans /= Math.pow(2, i - j);
		ans /= MathUtils.factorial(j);
		ans /= MathUtils.factorial(j - 1);

		return ans;
	}

	static double f_1(int x, int y) {
		double ans = 1.0;
		for (int i = 0; i < y; i++) {
			ans *= (x + i);
		}

		return ans;
	}

	static double f_2(int x, int y) {
		double ans = 1.0;
		for (int i = 0; i < y; i++) {
			ans *= (x - i);
		}

		return ans;
	}

	static long binomialInt(int n, int k) {
		if (k > n - k)
			k = n - k;

		long binom = 1;
		for (int i = 1; i <= k; i++)
			binom = binom * (n + 1 - i) / i;
		return binom;
	}

	static double gUp(int i, int j, double tau, double Ne) {
		double ans = 0.0;
		for (int k = j; k <= i; k++) {
			ans += (2 * k - 1) * Math.pow(-1, k - j) * f_1(j, k - 1) * f_2(i, k)
					* Math.exp(-(k * (k - 1) * tau * 0.5) / Ne) /
					(MathUtils.factorial(j) * MathUtils.factorial(k - j) * f_1(i, k));
		}

		return ans;
	}
}
		
		
		



	


