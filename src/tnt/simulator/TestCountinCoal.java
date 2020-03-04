package tnt.simulator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.util.Randomizer;

public class TestCountinCoal {

	public static void main(String args[]) {
		double bottleneckStrength = 1.0;
		double Ne = 1.0;
		int[] total_ways_to3 = new int[3]; // (4,1,1); (3,2,1); (2,2,2)
		int ways_to_321_in3 = 0;

		for (int j = 0; j < 100000000; j++) {
		int n=6;
		int[] ind = new int[n];
		for (int i=1; i<=6; i++) {
			ind[i-1]=i;
		}
		

		List<Integer[]> coal = new ArrayList<Integer[]>();
		
		for (int i = 1; i <= 6; i++) {
			coal.add(new Integer[] { i });
		}
		
		double duplicateTime = 0;
		double stopTime = 0 + bottleneckStrength;

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

				coal.add(concat(node1, node2));
				coal.remove(node1);
				coal.remove(node2);
				
				n = coal.size();
		}
			duplicateTime += deltaT;
		}
		
			List<Integer> lengths = new ArrayList<>();
		if (coal.size() == 3) {
			for (int i = 0; i < 3; i++) {
					lengths.add(coal.get(i).length);
				}
				if (lengths.contains(4))
					total_ways_to3[0] += 1;
				else if (lengths.contains(3)) {
					total_ways_to3[1] += 1;
					for (int i = 0; i < 3; i++) {
					List<Integer> list = Arrays.asList(coal.get(i));
					if (list.contains(1) && list.contains(2) && list.contains(3)) {
							for (int k = 0; k < 3 && k != i; k++) {
								List<Integer> list2 = Arrays.asList(coal.get(k));
								if (list2.contains(4) && list2.contains(5))
									ways_to_321_in3 += 1;
							}
					}
				}
				} else if (lengths.contains(2))
					total_ways_to3[2] += 1;
			}

		}
		System.out.println(Arrays.toString(total_ways_to3));
		System.out.println(ways_to_321_in3);
		System.out.println((double) ways_to_321_in3 / (double) total_ways_to3[1]);

	}

	public static <T> T[] concat(T[] first, T[] second) {
		T[] result = Arrays.copyOf(first, first.length + second.length);
		System.arraycopy(second, 0, result, first.length, second.length);
		return result;
	}
}
		
		
		



	


