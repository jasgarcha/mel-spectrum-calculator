package mel;

import discreteFourierTransform.DFT;

/**
 * MelSpectrum
 * @author Jasminder Garcha
 *
 */

public class MelSpectrum {
	/**
	 * The Mel Filter bank frequency table. The table determines which "bucket" to throw a frequency in.
	 * At index i, column 0 is the lower bound of the frequency range. Lower = CenterFrequency-(BandWidth/2). 
	 * At index i, column 1 is the upper bound of the frequency range. Upper = CenterFrequency+(BandWidth/2).
	 * From Professor Bon Sy's Mel-Spectrum-filter-band-frequency-range.xlsx	
	 */
	public final static double[][] MelScale = 
	{ 	//Lower, Upper. (Frequency (Hz)).
		{50, 150},
		{150, 250},
		{250, 350},
		{350, 450},
		{450, 550},
		{550, 650},
		{650, 750},
		{750, 850},
		{850, 950},
		{938, 1062},
		{1069, 1229},
		{1228, 1412},
		{1410.5, 1621.5},
		{1620, 1862},
		{1861, 2139},
		{2137, 2457}, 
		{2455.5, 2822.5},
		{2820, 3242},
		{3240, 3724},
		{3722, 4278}
	};		

	/**
	 * Filter function to determine if the frequency belongs in the "bucket" (frequency range).
	 * @param l The lth filter of the Mel Filter bank.
	 * @param frequency Continuous frequency in Hz.
	 * @return 1 if the frequency is within the specific "bucket", 0 if not due to the assumption of a square versus triangular filter. 
	 */
	public static int M(int l, double frequency) { 
		double lower = MelScale[l][0]; //index l of MelScale, column 0 is lower.
		double upper = MelScale[l][1]; //index l of MelScale, column 1 is upper.
		if(frequency >= lower && frequency <= upper) //Inclusive boundary. Check if the frequency is within range. 
			return 1;
		return 0;			
	}
	
	/**
	 * Computes the Mel Spectrum of a set of sampled data. 
	 * @param N The number of sampled data points. Only N/2 will be used because of symmetry. 
	 * @param samplingFrequency The sampling frequency in Hz.
	 * @param SOflReal Real term of Mel Spectrum S(l).
	 * @param SOflImaginary Imaginary term of Mel Spectrum S(l).
	 */
	public void computeMelSpectrum(int N, double samplingFrequency, double[] SOflReal, double[] SOflImaginary) {	
		int rangeOfk = (int)Math.ceil((double)N/2); //N/2. Round up for odd cases. 
		int rangeOfl = 20; //l=0,1,...,L-1. L is The total number of triangular Mel weighing filters (20).
		double[] continuousFrequency = new double[rangeOfk]; //Continuous frequency. 	
		double[] SOfkReal = new double[rangeOfk];
		double[] SOfkImaginary = new double[rangeOfk];
		for(int k = 0; k < continuousFrequency.length; k++) {
			continuousFrequency[k] = ((samplingFrequency/2)*k)/(rangeOfk); //Continuous frequency corresponds to a discrete frequency @ k. 
		}
		DFT dft = new DFT();
		dft.discreteFourierTransform(continuousFrequency, SOfkReal, SOfkImaginary); //Discrete Fourier transform of the continuous frequencies.  
		for(int l = 0; l < rangeOfl; l++) { //The range of l is from 0 to 19 (20-1).
			for(int k = 0; k < SOfkReal.length; k++) { //k: 0, 1,.., N/2.
				SOflReal[l] += SOfkReal[k]*M(l, continuousFrequency[k]);
				SOflImaginary[l] += SOfkImaginary[k]*M(l, continuousFrequency[k]);
			}
		}
	}
	
	/**
	 * Computes the Mel Spectrum of a set of sampled data.
	 * @param N The number of sampled data points. Only N/2 will be used because of symmetry. 
	 * @param samplingFrequency The sampling frequency in Hz.
	 * @return The Mel Spectrum is a 20x1 vector. The absolute value of each Mel Spectrum, |S(l)|, is stored in the rows.
	 */
	public double[][] computeMelSpectrum(int N, double samplingFrequency) {	
		int rangeOfk = (int)Math.ceil((double)N/2); //N/2. Round up for odd cases. 
		int rangeOfl = 20; //l=0,1,...,L-1. L is The total number of triangular Mel weighing filters (20).
		double[] continuousFrequency = new double[rangeOfk]; //Continuous frequency. 	
		double[] SOfkReal = new double[rangeOfk];
		double[] SOfkImaginary = new double[rangeOfk];
		double[] SOflReal = new double[20];
		double[] SOflImaginary = new double[20];
		double[][] MelSpectrum = new double[20][1]; //Mel Spectrum. 20x1 vector. 
		for(int k = 0; k < continuousFrequency.length; k++) 
			continuousFrequency[k] = ((samplingFrequency/2)*k)/(rangeOfk); //Continuous frequency corresponds to a discrete frequency @ k. 		
		DFT dft = new DFT();
		dft.discreteFourierTransform(continuousFrequency, SOfkReal, SOfkImaginary); //Discrete Fourier transform of the continuous frequencies. 
		for(int l = 0; l < rangeOfl; l++) { //The range of l is from 0 to 19 (20-1).
			for(int k = 0; k < SOfkReal.length; k++) { //k: 0, 1,.., N/2.
				SOflReal[l] += SOfkReal[k]*M(l, continuousFrequency[k]);
				SOflImaginary[l] += SOfkImaginary[k]*M(l, continuousFrequency[k]);
			}
			MelSpectrum[l][0] = Math.sqrt(SOflReal[l]*SOflReal[l]+SOflImaginary[l]*SOflImaginary[l]); //|S(l)|
		}
		return MelSpectrum;
	}
	
	/**
	 * Computes the Mel Spectrum of a set of sampled data and prints intermediate calculations. 
	 * @param N The number of sampled data points. Only N/2 will be used because of symmetry. 
	 * @param samplingFrequency The sampling frequency in Hz.
	 */
	public void printMelSpectrumWithSteps(int N, double samplingFrequency) {	
		int rangeOfk = (int)Math.ceil((double)N/2); //N/2. Round up for odd cases. 
		int rangeOfl = 20; //l=0,1,...,L-1. L is The total number of triangular Mel weighing filters (20).
		double[] continuousFrequency = new double[rangeOfk]; //Continuous frequency. 	
		double[] SOfkReal = new double[rangeOfk];
		double[] SOfkImaginary = new double[rangeOfk];
		double[] SOflReal = new double[20];
		double[] SOflImaginary = new double[20];
		String[][] steps = new String[20][2]; //Steps.
		for(int i = 0; i < steps.length; i++) { for(int j = 0; j < steps[i].length; j++) steps[i][j] = ""; } //Initialize null object to empty ("") string.
		for(int k = 0; k < continuousFrequency.length; k++) 
			continuousFrequency[k] = ((samplingFrequency/2)*k)/(rangeOfk); //Continuous frequency corresponds to a discrete frequency @ k. 		
		DFT dft = new DFT();
		dft.discreteFourierTransform(continuousFrequency, SOfkReal, SOfkImaginary); //Discrete Fourier transform of the continuous frequencies.  
		for(int l = 0; l < rangeOfl; l++) { //The range of l is from 0 to 19 (20-1).
			System.out.println("l="+l+":");
			System.out.println("S(l="+l+") = ");
			for(int k = 0; k < SOfkReal.length; k++) { //k: 0, 1,.., N/2.
				SOflReal[l] += SOfkReal[k]*M(l, continuousFrequency[k]);
				SOflImaginary[l] += SOfkImaginary[k]*M(l, continuousFrequency[k]);
				if(M(l, continuousFrequency[k]) == 1) { //M returns 1. The frequency belongs in this "bucket."
					steps[l][0] += "S(k="+k+")*1 + "; //Narrow down.
					steps[l][1] += "("+SOfkReal[k]+"+"+SOfkImaginary[k]+"j)*1 + ";
				}
				else {
					steps[l][0] += "S(k="+k+")*0 + ";
					steps[l][1] += "0 + ";
				}
			}
			//Intermediate calculations.
			for(int i = 0; i < steps[l].length; i++) { 
				steps[l][i] = steps[l][i].substring(0, steps[l][i].length()-3); //Truncate " + " on the last step.
				System.out.println(steps[l][i]);
			}
			System.out.println("S(l="+l+") = "+SOflReal[l]+"+"+SOflImaginary[l]+"j"); //S(l) = real term + imaginary term*j
			System.out.println();
		}
	}
}