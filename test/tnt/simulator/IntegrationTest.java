package tnt.simulator;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.junit.Assert;
import org.junit.Test;


public class IntegrationTest {

	@Test
	public void test() {
		double parentTime = 10;
		double lambda = 1.0;
		double mu = 0.5;
		double psi = 0.1;
		double currentTime = 0.1;
		double samplingExtantRate = 1.0;
		double delta = 1E-15;

		double u = 0.34;
		double c_1 = Math.abs(Math.sqrt(Math.pow((lambda - mu - psi), 2) + 4 * lambda * psi));
		double c_2 = -(lambda - mu - 2 * lambda * samplingExtantRate - psi) / c_1;


		BrentSolver solver = new BrentSolver(delta);

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
				return this.function(t) + Math.log(1 - u) - this.function(currentTime);
			}

			@Override
			public DerivativeStructure value(DerivativeStructure t) throws DimensionMismatchException {

				return this.ds(t).add(Math.log(1 - u)).subtract(this.ds(new DerivativeStructure(1, 1, currentTime)));

			}

		};
		double t_event = 1.8080588874991301;
		DerivativeStructure x = new DerivativeStructure(1, 1, t_event);

		
		Assert.assertEquals(f.value(x).getValue(), 0.0, delta);
		Assert.assertEquals(solver.solve(100000, f, currentTime, parentTime), t_event, delta);



	}

}
