/**
 * 
 */
package edu.cpp.cs.cs331;

import java.util.Random;
import java.util.Scanner;

/**
 * @author tarikrajper
 *
 */
public class Main {
	
	private static int[][] fir;	//First Input Matrix
	private static int[][] sec;	//Second Input Matrix
	private static int n;	//Size of Matrices
	
	private static Random rand;
	private static Random r;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		Scanner scan = new Scanner(System.in);
		
		System.out.println("Enter Matrix Size: ");	//User Input
		n = scan.nextInt();
		scan.close();
		
		System.out.println("-----------------------------------------");
		
		
		
		fir = new int[n][n];
		sec = new int[n][n];
		int[][] out = new int[n][n];		//Set user Input as dimensions
		
		rand = new Random(1000);
		r = new Random(1000);
		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				int x = rand.nextInt(10) + 1;
				int y = r.nextInt(10) + 1;
				fir[i][j] = x;
				sec[i][j] = y;
			}
		}
		long startTime = System.currentTimeMillis();						//Conducts Classical Multiplication and calculates running time
		classicalMultiplication(n, fir, sec, out);
		long endTime = System.currentTimeMillis();
		System.out.println("Classical Matrix Multiplication took: " + 
		(endTime - startTime) / 1 + " Milliseconds");
		printMatrix(out);
		System.out.println("---------------------------------------");
		
		long start = System.currentTimeMillis();							//Conducts Strassen's method and calculates running time
		out = strassen(fir, sec);
		long end = System.currentTimeMillis();
		System.out.println("Strassen's method took: " + 
		(end - start) / 1 + " Milliseconds");
		printMatrix(out);
		
	}
	
	/**
	 * Matrix multiplication (row by column)
	 * @param n Size of the matrix
	 * @param first	input
	 * @param second input
	 * @param result output
	 */
	public static void classicalMultiplication(int n, int[][] first, int[][] second, int[][] result) {
		for (int i = 0; i < n; i++) {	//Loop through rows
			for (int j = 0; j < n; j++) {	//Loop through Columns
				result[i][j] = 0;		//Initialize output
				
				for (int k = 0; k < n; k++) {
					result[i][j] += first[i][k] * second[k][j];	//Fill output with product of the two matrices
				}
			}
		}
	}
	
	/**
	 * 
	 * @param fir	Input
	 * @param sec	Second input
	 * @returns		Output matrix(product)
	 */
	public static int[][] strassen(int [][] fir, int[][] sec) {
		int n = fir.length;		//Length of the first Input
		
		int [][] result = new int[n][n];		//Sets dimensions of the output
		
		if ((n % 2 != 0) && (n != 1)) {		//Checks if the length is divisible by 2 and if the length is one
			int [][] a;
			int[][] b;		//Three temp matrices
			int [][] c;
			
			int len = n + 1;		//Length of the temp matrices (Add 1 for divisibility by 2)
			
			a = new int[len][len];
			b = new int[len][len];		//Initialize length of temp matrices
			c = new int[len][len];
			
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					a[i][j] = fir[i][j];		//Fill the temp matrices with values of the original inputs
					b[i][j] = sec[i][j];
				}
			}
			
			c = strassen(a, b);			//recursive call to implement strassens on the corrected sized matrices
			for (int i = 0; i < n; i++) 
				for (int j = 0; j < n; j++) 
					result[i][j] = c[i][j];		//Initialize result to the temp Matrix
				
			
			return result;
		}
		
		if (n == 1) {
			result[0][0] = fir[0][0] * sec[0][0];	//Normal Multiplication if n = 1
		} else {
			
			int[][] A11 = new int[n / 2][n / 2];		//initialize sub matrices for first input
			int[][] A12 = new int[n / 2][n / 2];
			int[][] A21 = new int[n / 2][n / 2];
			int[][] A22 = new int[n / 2][n / 2];
			
			int[][] B11 = new int[n / 2][n / 2];		//initialize sub matrices for second input
			int[][] B12 = new int[n / 2][n / 2];
			int[][] B21 = new int[n / 2][n / 2];
			int[][] B22 = new int[n / 2][n / 2];
			
			partition(fir, A11, 0, 0);				//Partitions the first Matrix into 4 sub matrices
			partition(fir, A12, 0, n / 2);
			partition(fir, A21, n / 2, 0);
			partition(fir, A22, n / 2, n / 2);
			
			partition(sec, B11, 0, 0);				//Partitions the second Matrix into 4 sub matrices
			partition(sec, B12, 0, n / 2);
			partition(sec, B21, n / 2, 0);
			partition(sec, B22, n / 2, n / 2);
			
			int[][] p = strassen(add(A11, A22), add(B11, B22)); //1			Performs operations on the sub matrices in accordance to strassens 
			int[][] q = strassen(add(A21, A22), B11); //2					method
			int[][] r = strassen(A11, subtract(B12, B22)); //3
			int[][] s = strassen(A22, subtract(B21, B11)); //4
			int[][] t = strassen(add(A11, A12), B22); //5
			int[][] u = strassen(subtract(A21, A11), add(B11, B12)); //6
			int[][] v = strassen(subtract(A12, A22), add(B21, B22)); //7
			
			int[][] C11 = add(subtract(add(p, s), t), v);					//Sets the sub matrices for the output to the values obtained from prev operations
			int[][] C12 = add(r, t);											//Part of Strassens method
			int[][] C21 = add(q, s);
			int[][] C22 = add(subtract(add(p, r), q), u);
			
			copyMatrix(C11, result, 0, 0);									//Copies the outputs sub matrices into the original output Matrix
			copyMatrix(C12, result, 0, n / 2);
			copyMatrix(C21, result, n / 2, 0);
			copyMatrix(C22, result, n / 2, n / 2);
			
		}
		return result;
	}
	
	public static int[][] add(int[][]A, int[][] B) {					//Adds two matrices together
		int n = A.length;
		
		int[][] result = new int[n][n];
		
		for (int i = 0; i < n; i++) 
			for (int j = 0; j < n; j++) 
				result[i][j] = A[i][j] + B[i][j];					//Sets sum equal to reult matrix
			
		
		return result;
	}
	
	public static int[][] subtract(int[][] A, int[][] B) {			//Subtracts two matrices together
		int n = A.length;
		
		int [][] result = new int[n][n];
		
		for (int i = 0; i < n; i++) 
			for (int j = 0; j < n; j++) 
				result[i][j] = A[i][j] - B[i][j];					//Sets difference into result matrix
			
		
		return result;
	}
	
	public static void partition(int[][] p, int[][] c, int k, int l) {		//Partitions the input matrices into smaller sub matrices
		for(int i0 = 0, i1= k; i0 < c.length; i0++, i1++) 
			for (int j0 = 0, j1 = l; j0 < c.length; j0++, j1++) 
				c[i0][j0] = p[i1][j1];	
	}
	
	public static void copyMatrix(int[][] c, int[][] p, int b, int c1) {		//Copies the data from the output's sub matrices into the original output Matrix
		for(int i1 = 0, i2 = b; i1 < c.length; i1++, i2++) 			
			for (int j1 = 0, j2 = c1; j1 < c.length; j1++, j2++) 
				p[i2][j2] = c[i1][j1];
			
		
	}
	
	
	public static void printMatrix(int[][] mat) {							//Prints Matrix
		for (int i = 0; i < mat.length; i++) {
			for (int j = 0; j < mat[i].length; j++) {
				System.out.print(mat[i][j] + " ");
			}
			System.out.println();
		}
	}
	

}
