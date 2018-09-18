
/**
 * Overall class that interpolates a cubic spline from a number of way points.
 * 
 * @author cooli
 *
 */
public class Spline {

	private float[] wayPoints;
	private int numWP;
	
	private float[] dxdtAtWP;
	private float[] dydtAtWP;
	
	private float[] ax;
	private float[] bx;
	private float[] ay;
	private float[] by;
	
	public Spline(float[] wayPoints) 
	{
		this.wayPoints = wayPoints;
		this.numWP = wayPoints.length / 3;
		
		this.dxdtAtWP = new float[numWP];
		this.dydtAtWP = new float[numWP];
		
		this.ax = new float[numWP];
		this.bx = new float[numWP];
		this.ay = new float[numWP];
		this.by = new float[numWP];
		
		create();
	}
	
	private int getIndex(float t)
	{
		for(int i=1; i<numWP; i++)
		{
			if(wayPoints[3*(i-1)+2] <= t && t < wayPoints[3*i+2])
			{
				return i;
			}
		}
		
		return -1;
	}

	public float getX(float t)
	{
		int i = getIndex(t);
		if(i == -1) return 0;
		
		float c = (t - wayPoints[3*(i-1)+2])
				/ (wayPoints[3*i+2]-wayPoints[3*(i-1)+2]);
		
		return (1-c)*wayPoints[3*(i-1)] + c*wayPoints[3*i] 
			 + c*(1-c)*(ax[i]*(1-c)+bx[i]*c);
	}
	
	public float getY(float t)
	{
		int i = getIndex(t);
		if(i == -1) return 0;
		
		float c = (t - wayPoints[3*(i-1)+2])
				/ (wayPoints[3*i+2]-wayPoints[3*(i-1)+2]);
		
		return (1-c)*wayPoints[3*(i-1)+1] + c*wayPoints[3*i+1] 
			 + c*(1-c)*(ay[i]*(1-c)+by[i]*c);
	}
	
	public float getVelX(float t)
	{
		int i = getIndex(t);
		if(i == -1) return 0;
		
		float c = (t - wayPoints[3*(i-1)+2])
				/ (wayPoints[3*i+2]-wayPoints[3*(i-1)+2]);
		
		float xi = wayPoints[3*i];
		float xim = wayPoints[3*(i-1)];
		float ti = wayPoints[3*i+2];
		float tim = wayPoints[3*(i-1)+2];
		
		return (xi - xim)/(ti - tim) + (1-2*c)*(ax[i]*(1-c)+bx[i]*c)/(ti - tim)
				+ c*(1-c)*(bx[i]-ax[i])/(ti - tim);
	}
	
	public float getVelY(float t)
	{
		int i = getIndex(t);
		if(i == -1) return 0;
		
		float c = (t - wayPoints[3*(i-1)+2])
				/ (wayPoints[3*i+2]-wayPoints[3*(i-1)+2]);
		
		float yi = wayPoints[3*i+1];
		float yim = wayPoints[3*(i-1)+1];
		float ti = wayPoints[3*i+2];
		float tim = wayPoints[3*(i-1)+2];
		
		return (yi - yim)/(ti - tim) + (1-2*c)*(ay[i]*(1-c)+by[i]*c)/(ti - tim)
				+ c*(1-c)*(by[i]-ay[i])/(ti - tim);
	}
	
	public float getAccX(float t)
	{
		int i = getIndex(t);
		if(i == -1) return 0;
		
		float c = (t - wayPoints[3*(i-1)+2])
				/ (wayPoints[3*i+2]-wayPoints[3*(i-1)+2]);
		
		float ti = wayPoints[3*i+2];
		float tim = wayPoints[3*(i-1)+2];
		
		return 2 * (bx[i] - 2*ax[i] + (ax[i]-bx[i])*3*c)
				 / ((ti-tim)*(ti-tim));
	}
	
	public float getAccY(float t)
	{
		int i = getIndex(t);
		if(i == -1) return 0;
		
		float c = (t - wayPoints[3*(i-1)+2])
				/ (wayPoints[3*i+2]-wayPoints[3*(i-1)+2]);
		
		float ti = wayPoints[3*i+2];
		float tim = wayPoints[3*(i-1)+2];
		
		return 2 * (by[i] - 2*ay[i] + (ay[i]-by[i])*3*c)
				 / ((ti-tim)*(ti-tim));
	}
	
	private void create()
	{	
		float[][] xMatrix = new float[numWP][numWP];
		float[]   xVector = new float[numWP];
		
		float[][] yMatrix = new float[numWP][numWP];
		float[]   yVector = new float[numWP];
		
		//inject n-1 mandatory equations into the matrix
		for(int i=1; i<numWP-1; i++)
		{
			float xi = wayPoints[3*i];
			float yi = wayPoints[3*i + 1];
			float ti = wayPoints[3*i + 2];
			
			float xim = wayPoints[3*(i-1)];
			float yim = wayPoints[3*(i-1) + 1];
			float tim = wayPoints[3*(i-1) + 2];
			
			float xip = wayPoints[3*(i+1)];
			float yip = wayPoints[3*(i+1) + 1];
			float tip = wayPoints[3*(i+1) + 2];
			
			xVector[i] = 3 * ( (xi - xim)/((ti - tim)*(ti - tim)) 
							 + (xip - xi)/((tip - ti)*(tip - ti)) );
			
			yVector[i] = 3 * ( (yi - yim)/((ti - tim)*(ti - tim)) 
					 		 + (yip - yi)/((tip - ti)*(tip - ti)) );
			
			xMatrix[i][i-1] = yMatrix[i][i-1] = 1/(yi - yim);
			xMatrix[i][i]   = yMatrix[i][i]   = 2 * ( 1/(yi - yim) + 1/(yip - yi) );
			xMatrix[i][i+1] = yMatrix[i][i+1] = 1/(yip - yi);
		}
		
		//boundary condition at beginning
		float x0 = wayPoints[0];
		float y0 = wayPoints[1];
		float t0 = wayPoints[2];
		
		float x1 = wayPoints[3];
		float y1 = wayPoints[4];
		float t1 = wayPoints[5];
		
		xMatrix[0][0] = 2 / (t1 - t0);
		xMatrix[0][1] = 1 / (t1 - t0);
		xVector[0] = 3 * (x1 - x0) / ( (t1-t0)*(t1-t0) );
		
		yMatrix[0][0] = 2 / (t1 - t0);
		yMatrix[0][1] = 1 / (t1 - t0);
		yVector[0] = 3 * (y1 - y0) / ( (t1-t0)*(t1-t0) );
		
		//boundary condition at end
		float xn = wayPoints[3*(numWP-1)];
		float yn = wayPoints[3*(numWP-1) + 1];
		float tn = wayPoints[3*(numWP-1) + 2];
		
		float xnm = wayPoints[3*(numWP-2)];
		float ynm = wayPoints[3*(numWP-2) + 1];
		float tnm = wayPoints[3*(numWP-2) + 2];
		
		xMatrix[numWP-1][numWP-2]   = 1 / (tn - tnm);
		xMatrix[numWP-1][numWP-1]   = 2 / (tn - tnm);
		xVector[numWP-1] = 3 * (xn - xnm) / ( (tn-tnm)*(tn-tnm) );
		
		yMatrix[numWP-1][numWP-2]   = 1 / (tn - tnm);
		yMatrix[numWP-1][numWP-1]   = 2 / (tn - tnm);
		yVector[numWP-1] = 3 * (yn - ynm) / ( (tn-tnm)*(tn-tnm) );
		
		dxdtAtWP = LinearSystem.solveLinearSystem(xMatrix, xVector);
		dydtAtWP = LinearSystem.solveLinearSystem(yMatrix, yVector);
		
		for(int i=1; i<numWP; i++)
		{
			ax[i] = dxdtAtWP[i-1]*(wayPoints[3*i+2]-wayPoints[3*(i-1)+2])
				  - (wayPoints[3*i]-wayPoints[3*(i-1)]);
			ay[i] = dydtAtWP[i-1]*(wayPoints[3*i+2]-wayPoints[3*(i-1)+2])
					  - (wayPoints[3*i+1]-wayPoints[3*(i-1)+1]);
			
			bx[i] = -dxdtAtWP[i]*(wayPoints[3*i+2]-wayPoints[3*(i-1)+1])
				  + (wayPoints[3*i]-wayPoints[3*(i-1)]);
			by[i] = -dydtAtWP[i]*(wayPoints[3*i+2]-wayPoints[3*(i-1)+1])
					  + (wayPoints[3*i+1]-wayPoints[3*(i-1)+1]);
		}
	}
	
}
