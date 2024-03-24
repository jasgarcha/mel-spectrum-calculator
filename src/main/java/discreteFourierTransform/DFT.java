package discreteFourierTransform;

import java.math.BigDecimal;
import java.math.RoundingMode;

/**
 * Discrete Fourier transform.
 * @author Jasminder Garcha
 */
public class DFT {
	/**
	 * Rounds a double to 3 decimal places for presentation, not intermediate calculation. 
	 * @param value The decimal number to be rounded up.
	 * @return The decimal number rounded up to three places. 
	 */
	public double round(double value) {
	    BigDecimal d = new BigDecimal(Double.toString(value));
	    d = d.setScale(3, RoundingMode.HALF_UP);
	    return d.doubleValue();
	}
	
	/**
	 * Computes the discrete Fourier transform of a continuous signal.
	 * @param sOfn s(n): a continuous signal over discrete time. The sampled data points. 
	 * @param SOfkReal The real terms of S(k): signal spectrum over discrete frequency domain. The length must equal the number of sampled data points in s(n).
	 * @param SOfkImaginary The imaginary terms of S(k): signal spectrum over discrete frequency domain. The length must equal the number of sampled data points in s(n).
	 */
	//Euler's Formula: e^(j*x)=cos(x)+j*sin(x), where j=sqrt(-1).
	public void discreteFourierTransform(double[] sOfn, double[] SOfkReal, double[] SOfkImaginary) {
		int N = sOfn.length; //Number of data points.
		for (int k = 0; k <= N-1; k++) { //Range of k: 0 to N-1.
			SOfkReal[k] = 0; //Real term.
			SOfkImaginary[k] = 0; //Imaginary term.
			for (int n = 0; n <= N-1; n++) { //n: 0 to N-1.
				double angle = -2*Math.PI*k*n/N;
				SOfkReal[k] += sOfn[n]*Math.cos(angle); 
				SOfkImaginary[k] += sOfn[n]*Math.sin(angle); 
			}
		}
	}	
	
	/**
	 * Computes the discrete Fourier transformation of a continuous signal with intermediate calculations. 
	 * @param sOfn s(n): a continuous signal over discrete time. The sampled data points. 
	 * @param SOfkReal The real terms of S(k): signal spectrum over discrete frequency domain. The length must equal the number of sampled data points in s(n).
	 * @param SOfkImaginary The imaginary terms of S(k): signal spectrum over discrete frequency domain. The length must equal the number of sampled data points in s(n).
	 * @param steps Intermediate calculations of computing the discrete Fourier transform.
	 */
	public void discreteFourierTransformWithSteps(double[] sOfn, double[] SOfkReal, double[] SOfkImaginary, String[][] steps) {
		int N = sOfn.length; //Number of data points.
		for(int i = 0; i < steps.length; i++) { for(int j = 0; j < steps[i].length; j++) steps[i][j] = ""; } //Initialize null object to empty ("") string.
		for (int k = 0; k <= N-1; k++) { //Range of k is 0 to N-1.
			SOfkReal[k] = 0; //Real term.
			SOfkImaginary[k] = 0; //Imaginary term.
			steps[k][0] += "k="+k+":\n"; //k.
			steps[k][0] += "S(k="+k+")=S("+k+"):\n"; //S(k) 
			for (int n = 0; n <= N-1; n++) { //n: 0 to N-1.
				double angle = -2*Math.PI*k*n/N;
				SOfkReal[k] += sOfn[n]*Math.cos(angle); 
				SOfkImaginary[k] += sOfn[n]*Math.sin(angle); 
				/*Steps: intermediate calculations*/
				steps[k][0] += "s("+n+")*e^(-2*PI*"+k+"*"+n+"/"+N+")*j + "; //DFT formula summation expansion. 
				steps[k][1] += "s("+n+")*cos(-2*PI*"+k+"*"+n+"/"+N+")+j*sin(-2*PI*"+k+"*"+n+"/"+N+") + ";  //Euler's formula substitution.
				steps[k][2] += "s("+n+")*cos(PI*"+(-2*k*n)+"/"+N+")+j*sin(PI*"+(-2*k*n)+"/"+N+") + "; //Grouping numerator, except PI. 
				if(-2*k*n == 0) //Check for division by zero in the numerator.
					steps[k][3]+= "s("+n+")*cos(0)+j*sin(0) + ";
				else 
					steps[k][3] += "s("+n+")*cos(PI*"+((double)-2*k*n/N)+")+j*sin(PI*"+((double)-2*k*n/N)+") + ";
				steps[k][4] += "s("+n+")*cos("+Math.round(Math.toDegrees(Math.PI*-2*k*n/N))+")+j*sin("+Math.round(Math.toDegrees(Math.PI*-2*k*n/N))+") + "; //Conversion from radians to degrees.
				//Round from double to int, if 0 or 1. Otherwise, maintain precision.
				if(round(Math.cos(angle)) == 0.0 || round(Math.cos(angle)) == -0.0 || round(Math.cos(angle)) == 1.0 || round(Math.cos(angle)) == -1.0) { 
					steps[k][5] += "s("+n+")*"+Math.round(round(Math.cos(angle)));
					steps[k][6] += sOfn[n]+"*["+Math.round(round(Math.cos(angle)));
					steps[k][7] += sOfn[n]*Math.round(round(Math.cos(angle)));
				}
				else { 
					steps[k][5] += "s("+n+")*"+round(Math.cos(angle));
					steps[k][6] += sOfn[n]+"*["+round(Math.cos(angle));
					steps[k][7] += sOfn[n]*round(Math.cos(angle));					
				}
				if(round(Math.sin(angle)) == 0.0 || round(Math.sin(angle)) == -0.0 || round(Math.sin(angle)) == 1.0 || round(Math.sin(angle)) == -1.0) {
					steps[k][5] += "+j*"+Math.round(round(Math.sin(angle)))+" + ";
					steps[k][6] += "+j*"+Math.round(round(Math.sin(angle)))+"] + ";
					steps[k][7] += "+j*"+sOfn[n]*Math.round(round(Math.sin(angle)))+" + ";
				}
				else {
					steps[k][5] += "+j*"+round(Math.sin(angle))+" + ";
					steps[k][6] += "+j*"+round(Math.sin(angle))+"] + ";
					steps[k][7] += "+j*"+sOfn[n]*round(Math.sin(angle))+" + ";
				}				
			}
			for(int i = 0; i < steps[k].length; i++) {
				steps[k][i] = steps[k][i].substring(0, steps[k][i].length()-3); //Truncate " + " on the last step.		
			}
			steps[k][7] += "\nS("+k+") = "+round(SOfkReal[k])+"+"+round(SOfkImaginary[k])+"j\n";			
		}		
	}
	
	/**
	 * Computes and prints the intermediate calculation of the discrete Fourier transformation. 
	 * @param sOfn s(n): a continuous signal over discrete time. The sampled data points. 
	 */
	public void printDiscreteFourierTransformWithSteps(double[] sOfn) {
		int N = sOfn.length; //Number of data points.
		String[] results = new String[N]; //k, S(k).
		String[][] steps = new String[N][8]; //Intermediate calculations.
		for(int i = 0; i < steps.length; i++) { for(int j = 0; j < steps[i].length; j++) steps[i][j] = ""; } //Initialize null object to empty ("") string.
		for (int k = 0; k <= N-1; k++) { //Range of k is 0 to N-1.
			double SOfkReal = 0; //Real term.
			double SOfkImaginary = 0; //Imaginary term.
			System.out.println("k="+k+":"); //k.
			System.out.println("S(k="+k+")=S("+k+"):"); //S(k) 
			for (int n = 0; n <= N-1; n++) { //n: 0 to N-1.
				double angle = -2*Math.PI*k*n/N;
				SOfkReal += sOfn[n]*Math.cos(angle);
				SOfkImaginary += sOfn[n]*Math.sin(angle); 
				/*Steps: intermediate calculations*/
				steps[k][0] += "s("+n+")*e^(-2*PI*"+k+"*"+n+"/"+N+")*j + "; //DFT formula summation expansion. 
				steps[k][1] += "s("+n+")*cos(-2*PI*"+k+"*"+n+"/"+N+")+j*sin(-2*PI*"+k+"*"+n+"/"+N+") + ";  //Euler's formula substitution.
				steps[k][2] += "s("+n+")*cos(PI*"+(-2*k*n)+"/"+N+")+j*sin(PI*"+(-2*k*n)+"/"+N+") + "; //Grouping numerator, except PI. 
				if(-2*k*n == 0) //Check for division by zero in the numerator.
					steps[k][3]+= "s("+n+")*cos(0)+j*sin(0) + ";
				else 
					steps[k][3] += "s("+n+")*cos(PI*"+((double)-2*k*n/N)+")+j*sin(PI*"+((double)-2*k*n/N)+") + ";
				steps[k][4] += "s("+n+")*cos("+Math.round(Math.toDegrees(Math.PI*-2*k*n/N))+")+j*sin("+Math.round(Math.toDegrees(Math.PI*-2*k*n/N))+") + "; //Conversion from radians to degrees.
				//Round from double to int, if 0 or 1. Otherwise, maintain precision.
				if(round(Math.cos(angle)) == 0.0 || round(Math.cos(angle)) == -0.0 || round(Math.cos(angle)) == 1.0 || round(Math.cos(angle)) == -1.0) {
					steps[k][5] += "s("+n+")*"+Math.round(round(Math.cos(angle)));
					steps[k][6] += sOfn[n]+"*["+Math.round(round(Math.cos(angle)));
					steps[k][7] += sOfn[n]*Math.round(round(Math.cos(angle)));
				}
				else { 
					steps[k][5] += "s("+n+")*"+round(Math.cos(angle));
					steps[k][6] += sOfn[n]+"*["+round(Math.cos(angle));
					steps[k][7] += sOfn[n]*round(Math.cos(angle));					
				}
				if(round(Math.sin(angle)) == 0.0 || round(Math.sin(angle)) == -0.0 || round(Math.sin(angle)) == 1.0 || round(Math.sin(angle)) == -1.0) {
					steps[k][5] += "+j*"+Math.round(round(Math.sin(angle)))+" + ";
					steps[k][6] += "+j*"+Math.round(round(Math.sin(angle)))+"] + ";
					steps[k][7] += "+j*"+sOfn[n]*Math.round(round(Math.sin(angle)))+" + ";
				}
				else {
					steps[k][5] += "+j*"+round(Math.sin(angle))+" + ";
					steps[k][6] += "+j*"+round(Math.sin(angle))+"] + ";
					steps[k][7] += "+j*"+sOfn[n]*round(Math.sin(angle))+" + ";
				}				
			}
			results[k] = sOfn[k]+", "+round(SOfkReal)+"+"+round(SOfkImaginary)+"j";
			//Output intermediate calculations for S(k).
			for(int i = 0; i < steps[k].length; i++) { 
				steps[k][i] = steps[k][i].substring(0, steps[k][i].length()-3); //Truncate " + " on the last step.
				System.out.println(steps[k][i]);
			}
			System.out.println("S("+k+") = "+round(SOfkReal)+"+"+round(SOfkImaginary)+"j"); //S(k)=real term+j*imaginary term.
			System.out.println();			
		}
		//Comma separated values.
		System.out.println("k, S(k)");
        //Print results.
        for (String result : results) System.out.println(result);
	}
		
	/**
	 * Computes the absolute value of a complex number.
	 * @param realTerm The real term of a complex number. a in a + bi.
	 * @param imaginaryTerm The imaginary term of a complex number. b in a + bi.
	 * @return The absolute value of a complex number: sqrt(realTerm^2+imaginaryTerm^2). |a + bi|=sqrt(a^2+b^2).
	 */
	public double absoluteValueOfComplexNumber(double realTerm, double imaginaryTerm) {
		return Math.sqrt(realTerm*realTerm+imaginaryTerm*imaginaryTerm);
	}
	
	/**
	 * Computes |S(k)|, absolute value of S(k).
	 * @param SOfkReal The real term of S(k).
	 * @param SOfkImaginary The imaginary term of S(k).
	 * @return The absolute value of S(k). |S(k)|=sqrt(SOfkReal^2+SOfkImaginary^2).
	 */
	public double[] absoluteValueOfDiscreteFourierTransform(double[] SOfkReal, double[] SOfkImaginary) {
		double[] abs = new double[SOfkReal.length];
		for(int i = 0; i < SOfkReal.length; i++)
			abs[i] = absoluteValueOfComplexNumber(SOfkReal[i], SOfkImaginary[i]);
		return abs; //|S(k)|
	}
	
	/**
	 * Computes and prints |S(k)|, the absolute value of S(k).
	 * @param SOfkReal The real term of S(k).
	 * @param SOfkImaginary SOfkImaginary The imaginary term of S(k).
	 */
	public void printAbsoluteValueOfDiscreteFourierTransform(double[] SOfkReal, double[] SOfkImaginary) {
		double[] absoluteValueSOfk = absoluteValueOfDiscreteFourierTransform(SOfkReal, SOfkImaginary);
		System.out.println("|S(k)|");
        for (double d : absoluteValueSOfk)
			System.out.println(round(d));
	}
}