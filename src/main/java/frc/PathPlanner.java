
/**
 * Path Planner class that uses an interpolating cubic spline to generate a path
 * between a number of way points. 
 * 
 * @author Ian Coolidge
 *
 */
public class PathPlanner {

	private static Spline spline;
	
	//waypoints are entered as (x1, y1, t1, x2, y2, t2,...), units are arbitrary.
	//TODO create a more dynamic system of entering waypoints
	private static float[] wayPoints = {0,    0,    0,
										1,    1,    1,
										2,    3,    2,
										0,    4,    3};
	
	public static void main(String[] args)
	{
		generateSpline();
	}
	
	private static void generateSpline() 
	{
		if(wayPoints.length % 3 == 0)
			spline = new Spline(wayPoints);
	}
	
}
