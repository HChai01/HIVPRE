package Feature_encoding;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.Writer;
import java.text.DecimalFormat;

/*
 * This Java code is designed to calculate the nucleotide composition in the coding sequence.
 * 
 * Feature list:
 * F1: Adenine_composition;
 * F2: Thymine_composition;
 * F3: Cytosine_composition;
 * F4: Guanine_composition;
 * F5: AT_composition;
 * F6: AC_composition;
 * F7: AG_composition;
 * F8: ApA_composition;
 * F9: ApT_composition;
 * F10: ApC_composition;
 * F11: ApG_composition;
 * F12: TpA_composition;
 * F13: TpT_composition;
 * F14: TpC_composition;
 * F15: TpG_composition;
 * F16: CpA_composition;
 * F17: CpT_composition;
 * F18: CpC_composition;
 * F19: CpG_composition;
 * F20: GpA_composition;
 * F21: GpT_composition;
 * F22: GpC_composition;
 * F23: GpG_composition.
 * 
 * Citation: Chai H, Gu Q, Hughes J, Robertson DL (2022) In silico prediction of HIV-1-host molecular interactions and their directionality. PLoS Comput Biol 18(2): e1009720. https://doi.org/10.1371/journal.pcbi.1009720.
 * 
 * */

public class Nucleotide_composition {
	// Path of the main folder
	static String topPath = "/Users/imac/Documents/HIVPRE/data/";
	
	static String CDSPath = topPath + "DNA_CDs_example.txt";
	static String outPath = topPath + "Nucleotide_composition.txt";
	
	static int readpointer;
	static int temporarypointer;
	static String strtxt = "";
	static String seq = "";
	
	static double basesum;
	static double[] ATCGNum = new double[5];
	
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
                		System.out.println("FASTA format error, Please check the headline of your CDS.");
                }else {
                		headinfo = lineTxt.trim().replace(">", "");;
                		strtxt = "S" + readpointer + "\t" + headinfo + "\t";
                		seq = "";
                		while((lineTxt = bufferedReader.readLine()) != null){
                			if(lineTxt.contains(">")) {
                				calculatenucleotide(seq);
                				calculatedinucleotide(seq);
                				write(outPath,strtxt+"\n");
                				readpointer = readpointer + 1;
                				headinfo = lineTxt.trim().replace(">", "");
                				strtxt = "S" + readpointer + "\t" + headinfo + "\t";
                        		seq = "";
                			}else {
                				seq = seq + lineTxt.trim();
                			}
                		}
                		calculatenucleotide(seq);
        				calculatedinucleotide(seq);
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
	
	public static void calculatenucleotide(String query){
		for(temporarypointer=0;temporarypointer<5;temporarypointer++) {
  			ATCGNum[temporarypointer] = 0;
  		}
		char a[] = query.toCharArray();
		for(temporarypointer=0;temporarypointer<a.length;temporarypointer++) {
			countATCG(a[temporarypointer]);
		}
		basesum = ATCGNum[1] + ATCGNum[2] + ATCGNum[3] + ATCGNum[4];
		strtxt = strtxt + decimalFormat.format(ATCGNum[1]/basesum) + "\t" +
				decimalFormat.format(ATCGNum[2]/basesum) + "\t" +
				decimalFormat.format(ATCGNum[3]/basesum) + "\t" +
				decimalFormat.format(ATCGNum[4]/basesum) + "\t" +
				decimalFormat.format((ATCGNum[1]+ATCGNum[2])/basesum) + "\t" +
				decimalFormat.format((ATCGNum[1]+ATCGNum[3])/basesum) + "\t" +
				decimalFormat.format((ATCGNum[1]+ATCGNum[4])/basesum) + "\t"; 
	}
	
	public static void calculatedinucleotide(String query){
		double count;
		char a[] = query.toCharArray();
		//ApA
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'A' && a[readpointer+1] == 'A') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		//ApT
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'A' && a[readpointer+1] == 'T') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		//ApC
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'A' && a[readpointer+1] == 'C') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		//ApG
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'A' && a[readpointer+1] == 'G') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		
		//TpA
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'T' && a[readpointer+1] == 'A') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		//TpT
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'T' && a[readpointer+1] == 'T') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		//TpC
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'T' && a[readpointer+1] == 'C') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		//TpG
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'T' && a[readpointer+1] == 'G') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		
		//CpA
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'C' && a[readpointer+1] == 'A') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		//CpT
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'C' && a[readpointer+1] == 'T') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		//CpC
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'C' && a[readpointer+1] == 'C') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		//count
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'C' && a[readpointer+1] == 'G') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		
		//GpA
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'G' && a[readpointer+1] == 'A') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		//GpT
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'G' && a[readpointer+1] == 'T') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		//GpC
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'G' && a[readpointer+1] == 'C') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "\t";
		//GpG
		count = 0;
		for(readpointer=0;readpointer<a.length-1;readpointer++) {
			if(a[readpointer] == 'G' && a[readpointer+1] == 'G') {
				count = count + 1;
			}
		}
		strtxt = strtxt + decimalFormat.format(100*count/(basesum-1)) + "";          
	}
	
	public static void countATCG(char base) {
		switch(base) {
			case	 'A':
				ATCGNum[1] = ATCGNum[1] + 1;
				break;
			case	 'T':
				ATCGNum[2] = ATCGNum[2] + 1;
				break;
			case	 'C':
				ATCGNum[3] = ATCGNum[3] + 1;
				break;
			case	 'G':
				ATCGNum[4] = ATCGNum[4] + 1;
				break;
			default:
				ATCGNum[0] = ATCGNum[0] + 1;
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
				"Adenine_composition" + "\t" + 
				"Thymine_composition" + "\t" +
				"Cytosine_composition" + "\t" +
				"Guanine_composition" + "\t" +
				"AT_composition" + "\t" +
				"AC_composition" + "\t" +
				"AG_composition" + "\t" +
				"ApA_composition" + "\t" +
				"ApT_composition" + "\t" +
				"ApC_composition" + "\t" +
				"ApG_composition" + "\t" +
				"TpA_composition" + "\t" +
				"TpT_composition" + "\t" +
				"TpC_composition" + "\t" +
				"TpG_composition" + "\t" +
				"CpA_composition" + "\t" +
				"CpT_composition" + "\t" +
				"CpC_composition" + "\t" +
				"CpG_composition" + "\t" +
				"GpA_composition" + "\t" +
				"GpT_composition" + "\t" +
				"GpC_composition" + "\t" +
				"GpG_composition" + "\n";
		write(outPath,info);
		readCDS(CDSPath);
	}
	
}
