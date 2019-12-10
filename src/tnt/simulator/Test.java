package tnt.simulator;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction;
import org.apache.commons.math3.analysis.solvers.NewtonRaphsonSolver;
import org.apache.commons.math3.exception.DimensionMismatchException;

public class Test {

	public static void main(String args[]) {
		double parentTime = 10;
		double lambda = 1.0;
		double mu = 0.5;
		double psi = 0.1;
		double currentTime = 0.1;

		double u = 0.34; // Randomizer.nextDouble();
		double c_1 = Math.abs(Math.sqrt(Math.pow((lambda - mu - psi), 2) + 4 * lambda * psi));
		double c_2 = -(lambda - mu - psi) / c_1;

		NewtonRaphsonSolver solver = new NewtonRaphsonSolver(1E-9);

		UnivariateDifferentiableFunction f = new UnivariateDifferentiableFunction() {

			private double function(double t) {
				return ((-2 * Math.log(Math.abs((c_2 - 1) * Math.exp(-c_1 * t) - c_2 - 1))
						+ (-c_1 * t) + (mu + psi + lambda) * t) / 2);

			}

			private DerivativeStructure ds(DerivativeStructure t) {
				return ((t.multiply(-c_1)).exp().multiply(c_2 - 1).subtract(c_2 + 1).abs().log().multiply(-2)
						.add(t.multiply(-c_1).exp().log()).add(t.multiply(lambda + mu + psi))).divide(2);
			}

			@Override
			public double value(double t) {
				System.out.println(this.function(currentTime));
				return this.function(t) + Math.log(1 - u) - this.function(currentTime);
			}

			@Override
			public DerivativeStructure value(DerivativeStructure t) throws DimensionMismatchException {

				return this.ds(t).add(Math.log(1 - u)).subtract(this.ds(new DerivativeStructure(1, 1, currentTime)));

			}

//			@Override
//			public double value(double t) {
//				return ((-2 * Math.log(Math.abs((c_2 - 1) * Math.exp(-c_1 * t) - c_2 - 1))
//						+ Math.log(Math.exp(-c_1 * t)) + (mu + psi + lambda) * t) / 2) + Math.log(1 - u) + Math.log(2);
//			}
//
//			@Override
//			public DerivativeStructure value(DerivativeStructure t) throws DimensionMismatchException {
//				return (((t.multiply(-c_1)).exp().multiply(c_2 - 1).subtract(c_2 + 1).abs().log().multiply(-2)
//						.add(t.multiply(-c_1).exp().log()).add(t.multiply(lambda + mu + psi))).divide(2))
//								.add(Math.log(1 - u)).add(Math.log(2));
//
//			}

		};
		DerivativeStructure x = new DerivativeStructure(1, 1, 0.5300252617550164);
		System.out.println(f.value(x).getValue());
		System.out.println(solver.solve(100000, f, currentTime, parentTime));
		
		


	}

}
