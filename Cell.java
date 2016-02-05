import java.io.*;
import java.util.*;
import java.lang.Math;
public class Cell {
    
    public static double[][] fillMatrix(String filename, int matrixSize, Point latticeA, Point latticeB, Point latticeC, int supercellA, int supercellB, int supercellC) throws IOException {
	
	double[][] atomMatrix = new double[matrixSize*supercellA*supercellB*supercellC][6];
	BufferedReader atomBuffer = new BufferedReader(new FileReader(filename));
	Scanner scan = new Scanner(atomBuffer);
	int countLine = 0;
	
	while(scan.hasNext()) {
	    String atomLine = scan.nextLine();
	    atomLine = atomLine.trim();
	    String[] atomLineList = atomLine.split(" ");
	    int atomNumber = Integer.parseInt(atomLineList[0]);
	    double xValue = Double.parseDouble(atomLineList[1]);
	    double yValue = Double.parseDouble(atomLineList[2]);
	    double zValue = Double.parseDouble(atomLineList[3]);
	    
	    for (int i=0; i<supercellA; i++) {
		for (int j=0; j<supercellB; j++) {
		    for (int k=0; k<supercellC; k++) {
			double atomX = xValue + i*latticeA.getX() + j*latticeB.getX() + k*latticeC.getX(); 
			double atomY = yValue + i*latticeA.getY() + j*latticeB.getY() + k*latticeC.getY(); 
			double atomZ = zValue + i*latticeA.getZ() + j*latticeB.getZ() + k*latticeC.getZ();
			atomMatrix[countLine][0] = atomX;
			atomMatrix[countLine][1] = atomY;
			atomMatrix[countLine][2] = atomZ;
			countLine++;
		    }
		}
	    }
	}
	return atomMatrix;
    }
    
    public static void writeMatrix(double[][] matrix) {
	int rows = matrix.length;
	for (int i=0; i<rows; i++) {
	    System.out.println((i+1) + " " + matrix[i][0] + " " + matrix[i][1] + " " + matrix[i][2] + " " + matrix[i][3] + " " + matrix[i][4] + " " + matrix[i][5]);
		
	}
    }
}