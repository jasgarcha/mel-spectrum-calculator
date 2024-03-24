import mel.MelSpectrum;

/**
 * Demonstrates the computeMelSpectrum and printMelSpectrumWithSteps methods of the MelSpectrum class using input of 10 data points and frequency of 2000 Hz.
 * @author Jasminder Garcha
 */
public class MelSpectrumExample {
	public static void main(String[] args) {
		int n = 10;
		int f = 2000;
		
		double[] SOflReal = new double[20];
		double[] SOflImaginary = new double[20];
		
		MelSpectrum melSpectrum = new MelSpectrum();

		melSpectrum.computeMelSpectrum(n, f, SOflReal, SOflImaginary);	
		for(int i = 0; i < SOflReal.length; i++)
			System.out.println(SOflReal[i]+"+"+SOflImaginary[i]+"j");
		System.out.println();
		
		double[][] spectrum = melSpectrum.computeMelSpectrum(n, f);
		for(int i = 0; i < spectrum.length; i++)
			for(int j = 0; j < spectrum[i].length; j++)
				System.out.println(spectrum[i][j]);
		System.out.println();
		
		melSpectrum.printMelSpectrumWithSteps(n, f);
	}
}