package Feature_encoding;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.Writer;
import java.text.DecimalFormat;

/*
 * This Java code is designed to calculate the amino acid composition in the protein sequence.
 * 
 * Feature list:
 * F1:	Alanine_composition
 * F2:	Arginine_composition
 * F3:	Asparagine_composition
 * F4:	Aspartic acid_composition
 * F5:	Cysteine_composition
 * F6:	Glutamic_acid_composition
 * F7:	Glutamine_composition
 * F8:	Glycine_composition
 * F9:	Histidine_composition
 * F10:	Isoleucine_composition
 * F11:	Leucine_composition
 * F12:	Lysine_composition
 * F13:	Methionine_composition
 * F14:	Phenylalanine_composition
 * F15:	Proline_composition
 * F16:	Serine_composition
 * F17:	Threonine_composition
 * F18:	Tryptophan_composition
 * F19:	Tyrosine_composition
 * F20:	Valine_composition
 * F21:	Aliphatic_composition(A/G/I/L/P/V)
 * F22:	Basic/Positive_charged_composition(R/H/K)
 * F23:	Sulfur(C/M)
 * F24:	Hydroxyl(S/T)
 * F25:	Acidic/Negative_charged(D/E)
 * F26:	Amide(N/Q)
 * F27:	Hydrophobic(A/C/I/L/M/F/W/V)
 * F28:	Neutral(G/H/P/S/T/Y)
 * F29:	Hydrophilic(R/N/D/Q/E/K)
 * F30:	Tiny(A/G/S)
 * F31:	Small(N/D/C/P/T)
 * F32:	Medium(Q/E/H/V)
 * F33:	Large(R/I/L/K/M)
 * F34:	Huge/Aromatic(F/W/Y)
 * F35:	Uncharged(A/N/C/Q/G/I/L/M/F/P/S/T/W/Y/V)
 * F36:	Polar(R/N/D/Q/E/H/K/S/T/Y)
 * F37:	Nonpolar(A/C/G/I/L/M/F/P/W/V)
 * 
 * Citation: Haiting Chai, Quan Gu, Joseph Hughes and David L. Robertson (202X), In silico prediction of HIV-1-host molecular interactions and their directionality, PLOS Computational Biology, XX(X): XXX-XXX, PMID: XXXXXXXX.
 * 
 * */

public class Amino_acid_composition {
	// Path of the main folder
	static String topPath = "/Users/imac/Documents/HIVPRE/data/";
	
	static String CDSPath = topPath + "Protein_example.txt";
	static String outPath = topPath + "Amino_acid_composition.txt";
	
	static int readpointer;
	static int temporarypointer;
	static String strtxt = "";
	static String seq = "";
	
	static double basesum;
	static double[] AminoNum = new double[21];
	
	static DecimalFormat decimalFormat = new DecimalFormat("0.0000");
	
	public static void readCDS(String Path){
		String headinfo;
		readpointer = 1;
		try {
			String encoding = "UTF-8";
            File listfile=new File(Path);
            if(listfile.isFile() && listfile.exists()){ 
            	InputStreamReader read = new InputStreamReader(
                new FileInputStream(listfile),encoding);                    
                BufferedReader bufferedReader = new BufferedReader(read);
                String lineTxt = "N/A";
                lineTxt = bufferedReader.readLine();
                if(!lineTxt.contains(">")) {
                		System.out.println("FASTA format error, Please check the headline of your protein sequence.");
                }else {
                		headinfo = lineTxt.trim().replace(">", "");;
                		strtxt = "S" + readpointer + "\t" + headinfo;
                		seq = "";
                		while((lineTxt = bufferedReader.readLine()) != null){
                			if(lineTxt.contains(">")) {
                				calculateaminoacid(seq);
                				calculateAAgroup(seq);
                				write(outPath,strtxt+"\n");
                				readpointer = readpointer + 1;
                				headinfo = lineTxt.trim().replace(">", "");
                				strtxt = "S" + readpointer + "\t" + headinfo;
                        		seq = "";
                			}else {
                				seq = seq + lineTxt.trim();
                			}
                		}
                		calculateaminoacid(seq);
                		calculateAAgroup(seq);
        				write(outPath,strtxt+"\n");
                }
                read.close();
            }else{
            	System.out.println("CDS File missing!");
            }
        } catch (Exception e) {
            System.out.println("CDS File Read error!");
            e.printStackTrace();
        }
    }
	
	public static void calculateaminoacid(String query){
		basesum = 0;
		for(temporarypointer=0;temporarypointer<21;temporarypointer++) {
			AminoNum[temporarypointer] = 0;
  		}
		char a[] = query.toCharArray();
		for(temporarypointer=0;temporarypointer<a.length;temporarypointer++) {
			countAmino(a[temporarypointer]);
		}
		for(temporarypointer=1;temporarypointer<21;temporarypointer++) {
			basesum = basesum + AminoNum[temporarypointer];
		}
		for(temporarypointer=1;temporarypointer<21;temporarypointer++) {
			strtxt = strtxt + "\t" + decimalFormat.format(AminoNum[temporarypointer]/basesum);
		}
	}
	
	public static void calculateAAgroup(String query){
//    	Aliphatic
//    	(A/G/I/L/P/V)
    		double Group1 = AminoNum[1] + AminoNum[8] + AminoNum[10] + AminoNum[11] + AminoNum[15] + AminoNum[20];
		strtxt = strtxt + "\t" + decimalFormat.format(Group1/basesum);
//    	Basic/Positive charged
//    	(R/H/K)
		double Group2 = AminoNum[2] + AminoNum[9] + AminoNum[12];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group2/basesum);
//    	Sulfur
//    	(C/M)
		double Group3 = AminoNum[5] + AminoNum[13];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group3/basesum);
//    	Hydroxyl
//    	(S/T)
		double Group4 = AminoNum[16] + AminoNum[17];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group4/basesum);
//    	Acidic/Negative charged 
//    	(D/E)
		double Group5 = AminoNum[4] + AminoNum[6];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group5/basesum);
//    	Amide
//    	(N/Q)
		double Group6 = AminoNum[3] + AminoNum[7];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group6/basesum);
//    	Hydrophobic
//    	(A/C/I/L/M/F/W/V)
		double Group7 = AminoNum[1] + AminoNum[5] + AminoNum[10] + AminoNum[11] + AminoNum[13] + AminoNum[14] + AminoNum[18] + AminoNum[20];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group7/basesum);
//    	Neutral
//    	(G/H/P/S/T/Y)
		double Group8 = AminoNum[8] + AminoNum[9] + AminoNum[15] + AminoNum[16] + AminoNum[17] + AminoNum[19];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group8/basesum);
//    	Hydrophilic
//    	(R/N/D/Q/E/K)
		double Group9 = AminoNum[2] + AminoNum[3] + AminoNum[4] + AminoNum[7] + AminoNum[6] + AminoNum[12];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group9/basesum);
//    	Tiny
//    	(A/G/S)
		double Group10 = AminoNum[1] + AminoNum[8] + AminoNum[16];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group10/basesum);
//    	Small
//    	(N/D/C/P/T)
		double Group11 = AminoNum[3] + AminoNum[4] + AminoNum[5] + AminoNum[15] + AminoNum[17];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group11/basesum);
//    	Medium
//    	(Q/E/H/V)
		double Group12 = AminoNum[7] + AminoNum[6] + AminoNum[9] + AminoNum[20];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group12/basesum);
//    	Large
//    	(R/I/L/K/M)
		double Group13 = AminoNum[2] + AminoNum[10] + AminoNum[11] + AminoNum[12] + AminoNum[13];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group13/basesum);
//    	Huge/Aromatic
//    	(F/W/Y)
		double Group14 = AminoNum[14] + AminoNum[18] + AminoNum[19];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group14/basesum);
//    	Uncharged
//    	(A/N/C/Q/G/I/L/M/F/P/S/T/W/Y/V)
		double Group15 = AminoNum[1] + AminoNum[3] + AminoNum[5] + AminoNum[7] + AminoNum[8] + AminoNum[10] + AminoNum[11] + AminoNum[13] + AminoNum[14] + AminoNum[15] + AminoNum[16] + AminoNum[17] + AminoNum[18] + AminoNum[19] + AminoNum[20];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group15/basesum);
//    	Polar
//    	(R/N/D/Q/E/H/K/S/T/Y)
		double Group16 = AminoNum[2] + AminoNum[3] + AminoNum[4] + AminoNum[7] + AminoNum[6] + AminoNum[9] + AminoNum[12] + AminoNum[16] + AminoNum[17] + AminoNum[19];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group16/basesum);
//    	Nonpolar
//    	(A/C/G/I/L/M/F/P/W/V)
		double Group17 = AminoNum[1] + AminoNum[5] + AminoNum[8] + AminoNum[10] + AminoNum[11] + AminoNum[13] + AminoNum[14] + AminoNum[15] + AminoNum[18] + AminoNum[20];
		strtxt = strtxt +  "\t" + decimalFormat.format(Group17/basesum);		
	}
	
	public static void countAmino(char amino) {
		switch(amino) {
			case	 'A':
				AminoNum[1] = AminoNum[1] + 1;
				break;
			case	 'R':
				AminoNum[2] = AminoNum[2] + 1;
				break;
			case	 'N':
				AminoNum[3] = AminoNum[3] + 1;
				break;
			case	 'D':
				AminoNum[4] = AminoNum[4] + 1;
				break;
			case	 'C':
				AminoNum[5] = AminoNum[5] + 1;
				break;
			case	 'E':
				AminoNum[6] = AminoNum[6] + 1;
				break;
			case	 'Q':
				AminoNum[7] = AminoNum[7] + 1;
				break;
			case	 'G':
				AminoNum[8] = AminoNum[8] + 1;
				break;
			case	 'H':
				AminoNum[9] = AminoNum[9] + 1;
				break;
			case	 'I':
				AminoNum[10] = AminoNum[10] + 1;
				break;
			case	 'L':
				AminoNum[11] = AminoNum[11] + 1;
				break;
			case	 'K':
				AminoNum[12] = AminoNum[12] + 1;
				break;
			case	 'M':
				AminoNum[13] = AminoNum[13] + 1;
				break;
			case	 'F':
				AminoNum[14] = AminoNum[14] + 1;
				break;
			case	 'P':
				AminoNum[15] = AminoNum[15] + 1;
				break;
			case	 'S':
				AminoNum[16] = AminoNum[16] + 1;
				break;
			case	 'T':
				AminoNum[17] = AminoNum[17] + 1;
				break;
			case	 'W':
				AminoNum[18] = AminoNum[18] + 1;
				break;
			case	 'Y':
				AminoNum[19] = AminoNum[19] + 1;
				break;
			case	 'V':
				AminoNum[20] = AminoNum[20] + 1;
				break;
			default:
				AminoNum[0] = AminoNum[0] + 1;
				break;
		}
	}
	
	public static void write(String Path, String textstr){
		try {
			File textfile = new File(Path);	
			Writer output = null;
			output = new FileWriter(textfile,true);
			output.write(textstr);
			output.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	
	
	public static void main(String argv[]){
		String info;
		info = "Index" + "\t" + "Headinfo" + "\t" + 
				"Alanine_composition" + "\t" + 
				"Arginine_composition" + "\t" +
				"Asparagine_composition" + "\t" +
				"Aspartic acid_composition" + "\t" +
				"Cysteine_composition" + "\t" +
				"Glutamic_acid_composition" + "\t" +
				"Glutamine_composition" + "\t" +
				"Glycine_composition" + "\t" +
				"Histidine_composition" + "\t" +
				"Isoleucine_composition" + "\t" +
				"Leucine_composition" + "\t" +
				"Lysine_composition" + "\t" +
				"Methionine_composition" + "\t" +
				"Phenylalanine_composition" + "\t" +
				"Proline_composition" + "\t" +
				"Serine_composition" + "\t" +
				"Threonine_composition" + "\t" +
				"Tryptophan_composition" + "\t" +
				"Tyrosine_composition" + "\t" +
				"Valine_composition" + "\t" +
				"Aliphatic_composition(A/G/I/L/P/V)" + "\t" +
				"Basic/Positive_charged_composition(R/H/K)" + "\t" +
				"Sulfur(C/M)" + "\t" +
				"Hydroxyl(S/T)" + "\t" +
				"Acidic/Negative_charged(D/E)" + "\t" +
				"Amide(N/Q)" + "\t" +
				"Hydrophobic(A/C/I/L/M/F/W/V)" + "\t" +
				"Neutral(G/H/P/S/T/Y)" + "\t" +
				"Hydrophilic(R/N/D/Q/E/K)" + "\t" +
				"Tiny(A/G/S)" + "\t" +
				"Small(N/D/C/P/T)" + "\t" +
				"Medium(Q/E/H/V)" + "\t" +
				"Large(R/I/L/K/M)" + "\t" +
				"Huge/Aromatic(F/W/Y)" + "\t" +
				"Uncharged(A/N/C/Q/G/I/L/M/F/P/S/T/W/Y/V)" + "\t" +
				"Polar(R/N/D/Q/E/H/K/S/T/Y)" + "\t" +
				"Nonpolar(A/C/G/I/L/M/F/P/W/V)" + "\n";
		write(outPath,info);
		readCDS(CDSPath);
	}
	
}
