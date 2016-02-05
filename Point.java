import java.lang.Math;
public class Point {

    private double xLoc, yLoc, zLoc;

    public Point(double x, double y, double z) {
	xLoc = x;
	yLoc = y;
	zLoc = z;
    }
    
    public Point(Point p1) {
	xLoc = p1.getX();
	yLoc = p1.getY();
	zLoc = p1.getZ();
    }

    public Point() {
	xLoc = 0.0;
	yLoc = 0.0;
	zLoc = 0.0;
    }
    
    public void setX(double x) {
	xLoc = x;
    }

    public void setY(double y) {
	yLoc = y;
    }

    public void setZ(double z) {
	zLoc = z;
    }

    public double getX() {
	return xLoc;
    }
    public double getY() {
	return yLoc;
    }
    public double getZ() {
	return zLoc;
    }

    public double fromOrigin(){
	return Math.sqrt(xLoc*xLoc + yLoc*yLoc + zLoc*zLoc);
    }

    public Point multiply(double n) {
	double xNew = xLoc*n;
	double yNew = yLoc*n;
	double zNew = zLoc*n;
	return new Point(xNew, yNew, zNew);
    }


    public Point multiply(int n) {
	double xNew = xLoc*n;
	double yNew = yLoc*n;
	double zNew = zLoc*n;
	return new Point(xNew, yNew, zNew);
    }
    
    public static double magnitude(Point p1) {
	double x = p1.getX(); 
	double y = p1.getY(); 
	double z = p1.getZ();
	return Math.sqrt(x*x + y*y + z*z);
    }

    public static double separation(Point p1, Point p2) {
	double xDelta = p1.getX() - p2.getX();
	double yDelta = p1.getY() - p2.getY();
	double zDelta = p1.getZ() - p2.getZ();
	return Math.sqrt(xDelta*xDelta + yDelta*yDelta + zDelta*zDelta);
    }

    public static double separation2(Point p1, Point p2) {
	Point p3 = subtract(p1, p2);
	return magnitude(p3);
    }

    public static Point add(Point P1, Point P2) {
	double xNew = P1.getX() + P2.getX();
	double yNew = P1.getY() + P2.getY();
	double zNew = P1.getZ() + P2.getZ();
	return new Point(xNew, yNew, zNew);
    }

    public static Point subtract(Point P1, Point P2) {
	double xNew = P1.getX() - P2.getX();
	double yNew = P1.getY() - P2.getY();
	double zNew = P1.getZ() - P2.getZ();
	return new Point(xNew, yNew, zNew);
    }
    
    public static Point add(Point P1, Point P2, Point P3) {
	double xNew = P1.getX() + P2.getX() + P3.getX();
	double yNew = P1.getY() + P2.getY() + P3.getY();
	double zNew = P1.getZ() + P2.getZ() + P3.getZ();
	return new Point(xNew, yNew, zNew);
    }
    
    public static double minImage(Point p1, Point p2, Point latticeA, Point latticeB, Point latticeC) {
        Point r = subtract(p1, p2);
        double minImageDistance = magnitude(r);
	
	for (int i=-1; i<2; i++) {
	    for (int j=-1; j<2; j++) {
		for (int k=-1; k<2; k++) {
		    double imageX = r.getX() + i*latticeA.getX() + j*latticeB.getX() + k*latticeC.getX(); 
		    double imageY = r.getY() + i*latticeA.getY() + j*latticeB.getY() + k*latticeC.getY(); 
		    double imageZ = r.getZ() + i*latticeA.getZ() + j*latticeB.getZ() + k*latticeC.getZ();
		    
		    Point imageVector = new Point(imageX, imageY, imageZ);
		    double imageDistance = magnitude(imageVector);
		    
		    if (imageDistance < minImageDistance) {
			minImageDistance = imageDistance;
		    }
		}
	    }
	}
	return minImageDistance;
    }
    
    public static Point minImageVector(Point p1, Point p2, Point boxA, Point boxB, Point boxC) {
	Point r = subtract(p1, p2);
	double minImageDistance = magnitude(r);
	Point minImageVector = new Point(r);
	int degen = 0;
	for (int i=-1; i<2; i++) {
	    for (int j=-1; j<2; j++) {
		for (int k=-1; k<2; k++) {
		    double imageX = r.getX() + i*boxA.getX() + j*boxB.getX() + k*boxC.getX(); 
		    double imageY = r.getY() + i*boxA.getY() + j*boxB.getY() + k*boxC.getY(); 
		    double imageZ = r.getZ() + i*boxA.getZ() + j*boxB.getZ() + k*boxC.getZ();
		    
		    Point imageVector = new Point(imageX, imageY, imageZ);
		    double imageDistance = magnitude(imageVector);
		    
		    if (imageDistance < minImageDistance) {
			minImageDistance = imageDistance;
			minImageVector.setX(imageX); 
			minImageVector.setY(imageY);
			minImageVector.setZ(imageZ); 
		    }
		    
		    if (imageDistance==minImageDistance) {
			degen++;
			//	System.out.println("IMAGE DEGENERACY! = " + degen);
		    }
		}
	    }
	}
	return minImageVector;
    }
    
    public String toString() {
	return String.format("(%g,%g,%g)",xLoc,yLoc,zLoc);
    }

    public void write() {
	System.out.println(String.format("x: %g, y: %g,z: %g, mag: %g",xLoc,yLoc,zLoc,Math.sqrt(xLoc*xLoc + yLoc*yLoc + zLoc*zLoc)));
    }
    
}