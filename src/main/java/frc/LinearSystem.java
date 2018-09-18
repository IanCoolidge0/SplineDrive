package frc;

/**
 * Utility class based on 2d arrays to solve nxn linear systems.
 * 
 * @author cooli
 *
 */
public class LinearSystem {

	/**
	 * Solve an nxn linear system using Cramer's rule.
	 * 
	 * @param A the matrix in the equation Ax = b
	 * @param b the resultant vector in the equation Ax = b
	 * @return the resultant vector if the system is nondegenerate, null otherwise 
	 */
	public static float[] solveLinearSystem(float[][] A, float[] b)
	{
		if(!(A.length == A[0].length && A.length == b.length))
			return null;
		
		float[] out = new float[A.length];
		float[] old = new float[A.length];
		float detA = determinant(A);
		
		if(detA == 0)
			return null;
		
		for(int i=0; i<A.length; i++)
		{
			for(int j=0; j<A.length; j++) {
				old[j] = A[j][i];
				A[j][i] = b[j];
			}
			out[i] = determinant(A)/detA;
			for(int j=0; j<A.length; j++)
				A[j][i] = old[j];
		}
		
		return out;
	}
	
	/**
	 * Recursively calcualate the determinant of an nxn matrix A
	 * 
	 * @param A the matrix to calculate the determinant of
	 * @return the determinant
	 */
	public static float determinant(float[][] A)
	{
		if(A.length != A[0].length) {
			return 0;
		}
				
		if(A.length == 1) {
			return A[0][0];
		}
		
		float det = 0;
		
		for(int i=0; i<A.length; i++)
		{
			float[][] minor = new float[A.length-1][];
			
			for(int j1=0; j1<A.length-1; j1++)
				minor[j1] = new float[A.length-1];

			int rowIndex = -1;
			boolean skipped = false;
			for(int j=0; j<A.length; j++)
			{
				rowIndex++;
				if(i == j) {
					skipped = true;
					continue;
				}
				
				for(int k=0; k<A.length-1; k++)
				{
					minor[(skipped?j-1:j)][k] = A[k+1][rowIndex];
				}
			}
			
			det += (1-2*(i%2)) * A[0][i] * determinant(minor);
		}
		
		return det;
	}
	
}
