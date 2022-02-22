package Feature_encoding;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.Writer;
import java.text.DecimalFormat;

/*
 * This Java code is designed to calculate the codon usage in the mRNA.
 * 
 * Feature list:
 * F1: UUU_(F)
 * F2: UUC_(F)
 * F3: UUA_(L)
 * F4: UUG_(L)
 * F5: CUU_(L)
 * F6: CUC_(L)
 * F7: CUA_(L)
 * F8: CUG_(L)
 * F9: AUU_(I)
 * F10: AUC_(I)
 * F11: AUA_(I)
 * F12: AUG_(M)
 * F13: GUU_(V)
 * F14: GUC_(V)
 * F15: GUA_(V)
 * F16: GUG_(V)
 * F17: UCU_(S)
 * F18: UCC_(S)
 * F19: UCA_(S)
 * F20: UCG_(S)
 * F21: AGU_(S)
 * F22: AGC_(S)
 * F23: CCU_(P)
 * F24: CCC_(P)
 * F25: CCA_(P)
 * F26: CCG_(P)
 * F27: ACU_(T)
 * F28: ACC_(T)
 * F29: ACA_(T)
 * F30: ACG_(T)
 * F31: UAU_(Y)
 * F32: UAC_(Y)
 * F33: GCU_(A)
 * F34: GCC_(A)
 * F35: GCA_(A)
 * F36: GCG_(A)
 * F37: CAU_(H)
 * F38: CAC_(H)
 * F39: CAA_(Q)
 * F40: CAG_(Q)
 * F41: AAU_(N)
 * F42: AAC_(N)
 * F43: AAA_(K)
 * F44: AAG_(K)
 * F45: GAU_(D)
 * F46: GAC_(D)
 * F47: GAA_(E)
 * F48: GAG_(E)
 * F49: UGU_(C)
 * F50: UGC_(C)
 * F51: UGG_(W)
 * F52: CGU_(R)
 * F53: CGC_(R)
 * F54: CGA_(R)
 * F55: CGG_(R)
 * F56: AGA_(R)
 * F57: AGG_(R)
 * F58: GGU_(G)
 * F59: GGC_(G)
 * F60: GGA_(G)
 * F61: GGG_(G)
 * F62: UAA_(N/A)
 * F63: UAG_(N/A)
 * F64: UGA_(N/A)
 * 
 * Citation: Chai H, Gu Q, Hughes J, Robertson DL (2022) In silico prediction of HIV-1-host molecular interactions and their directionality. PLoS Comput Biol 18(2): e1009720. https://doi.org/10.1371/journal.pcbi.1009720.
 * 
 * */

public class Codon_usage_mRNA {
	// Path of the main folder
	static String topPath = "/Users/imac/Documents/HIVPRE/data/";
	
	static String mRNAPath = topPath + "mRNA_example.txt";
	static String outPath = topPath + "Condon_usage(ver.mRNA).txt";
	
	static int readpointer;
	static int temporarypointer;
	static String strtxt = "";
	static String seq = "";
	static int[] codonINFO = new int[65];
	
	static DecimalFormat decimalFormat = new DecimalFormat("0.0000");
	
	public static void readmRNA(String Path){
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
                		System.out.println("Please check the headline of your mRNA.");
                }else {
                		headinfo = lineTxt.trim().replace(">", "");;
                		strtxt = "S" + readpointer + "\t" + headinfo;
                		seq = "";
                		while((lineTxt = bufferedReader.readLine()) != null){
                			if(lineTxt.contains(">")) {
                				write(outPath,strtxt);
                				calculatecodon(seq,headinfo);              				
                				readpointer = readpointer + 1;
                				headinfo = lineTxt.trim().replace(">", "");
                				strtxt = "S" + readpointer + "\t" + headinfo;
                        		seq = "";
                			}else {
                				seq = seq + lineTxt.trim();
                			}
                		}
                		write(outPath,strtxt);
                		calculatecodon(seq,headinfo);
                }
                read.close();
            }else{
            	System.out.println("mRNA File missing!");
            }
        } catch (Exception e) {
            System.out.println("mRNA File Read error!");
            e.printStackTrace();
        }
    }
	
	public static void calculatecodon(String query,String id){
		String stopcodon;
		String startcodon;
		String clue = "";
		int key = 0;
		for(temporarypointer=0;temporarypointer<65;temporarypointer++) {
			codonINFO[temporarypointer] = 0;
  		}
		char a[] = query.toCharArray();
		if(a.length%3!=0 || a.length <6) {
			System.out.println("mRNA is incompleted, please check the content in: " + id);
			write(outPath,"\tmRNA incompleted\n");
		}else {
			stopcodon = a[a.length-3] + "" + a[a.length-2] + "" + a[a.length-1] + "";
			if(!(stopcodon.equals("UAA")||stopcodon.equals("UAG")||stopcodon.equals("UGA"))) {
				System.out.println("Stop codon is missing, please check the content in: " + id);
				write(outPath,"\tStop codon missing\n");
			}else {
				startcodon = a[0] + "" + a[1] + "" + a[2] + "";
				if(!startcodon.equals("AUG")) {
					System.out.println("Start codon is abnormal, please check the content in: " + id);
					write(outPath,"\tStart codon abnormal");
				}
				for(temporarypointer=0;temporarypointer<a.length/3;temporarypointer++) {
	        			clue = a[temporarypointer*3] + "" + a[temporarypointer*3+1] + "" + a[temporarypointer*3+2] + "";
	        			countcodon(clue);
	        			if(temporarypointer < (a.length/3)-1 && (clue.equals("UAA")||clue.equals("UAG")||clue.equals("UGA"))) {
	        				System.out.println("Abnormal stop codons '" + clue + "' is detected, please check the content in: " + id);
	        				key = 1;
	        			}
        			}
				if(key == 1) {
					write(outPath,"\tMultiple stop codons detected");
				}else {
					write(outPath,"\tNormal");
				}
				
				for(temporarypointer=1;temporarypointer<65;temporarypointer++) {
	        			if(codonINFO[temporarypointer]==0) {
	        				write(outPath,"\t0.0000");
	        			}else {          			
	                		if(temporarypointer>=1 && temporarypointer<=2) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[1]+codonINFO[2])));
	                		}
	                		if(temporarypointer>=3 && temporarypointer<=8) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[3]+codonINFO[4]+codonINFO[5]+codonINFO[6]+codonINFO[7]+codonINFO[8])));
	                		}
	                		if(temporarypointer>=9 && temporarypointer<=11) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[9]+codonINFO[10]+codonINFO[11])));
	                		}
	                		if(temporarypointer==12) {
	                			write(outPath,"\t1.0000");
	                		}
	                		if(temporarypointer>=13 && temporarypointer<=16) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[13]+codonINFO[14]+codonINFO[15]+codonINFO[16])));
	                		}
	                		if(temporarypointer>=17 && temporarypointer<=22) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[17]+codonINFO[18]+codonINFO[19]+codonINFO[20]+codonINFO[21]+codonINFO[22])));
	                		}
	                		if(temporarypointer>=23 && temporarypointer<=26) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[23]+codonINFO[24]+codonINFO[25]+codonINFO[26])));
	                		}
	                		if(temporarypointer>=27 && temporarypointer<=30) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[27]+codonINFO[28]+codonINFO[29]+codonINFO[30])));
	                		}
	                		if(temporarypointer>=31 && temporarypointer<=32) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[31]+codonINFO[32])));
	                		}
	                		if(temporarypointer>=33 && temporarypointer<=36) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[33]+codonINFO[34]+codonINFO[35]+codonINFO[36])));
	                		}
	                		if(temporarypointer>=37 && temporarypointer<=38) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[37]+codonINFO[38])));
	                		}
	                		if(temporarypointer>=39 && temporarypointer<=40) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[39]+codonINFO[40])));
	                		}
	                		if(temporarypointer>=41 && temporarypointer<=42) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[41]+codonINFO[42])));
	                		}
	                		if(temporarypointer>=43 && temporarypointer<=44) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[43]+codonINFO[44])));
	                		}
	                		if(temporarypointer>=45 && temporarypointer<=46) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[45]+codonINFO[46])));
	                		}
	                		if(temporarypointer>=47 && temporarypointer<=48) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[47]+codonINFO[48])));
	                		}
	                		if(temporarypointer>=49 && temporarypointer<=50) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[49]+codonINFO[50])));
	                		}
	                		if(temporarypointer==51) {
	                			write(outPath,"\t1.0000");
	                		}
	                		if(temporarypointer>=52 && temporarypointer<=57) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[52]+codonINFO[53]+codonINFO[54]+codonINFO[55]+codonINFO[56]+codonINFO[57])));
	                		}
	                		if(temporarypointer>=58 && temporarypointer<=61) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[58]+codonINFO[59]+codonINFO[60]+codonINFO[61])));
	                		}
	                		if(temporarypointer>=62 && temporarypointer<=64) {
	                			write(outPath,"\t" +decimalFormat.format((float)codonINFO[temporarypointer]/(codonINFO[62]+codonINFO[63]+codonINFO[64])));
	                		}
	        			}
				}
				write(outPath,"\n");
			}
		}		
	}
	
	public static void countcodon(String codon) {
		switch(codon) {
		case	 "UUU":
			codonINFO[1] = codonINFO[1] + 1;
			break;
		case	 "UUC":
			codonINFO[2] = codonINFO[2] + 1;
			break;
		case	 "UUA":
			codonINFO[3] = codonINFO[3] + 1;
			break;
		case	 "UUG":
			codonINFO[4] = codonINFO[4] + 1;
			break;
		case	 "CUU":
			codonINFO[5] = codonINFO[5] + 1;
			break;
		case	 "CUC":
			codonINFO[6] = codonINFO[6] + 1;
			break;
		case	 "CUA":
			codonINFO[7] = codonINFO[7] + 1;
			break;
		case	 "CUG":
			codonINFO[8] = codonINFO[8] + 1;
			break;
		case	 "AUU":
			codonINFO[9] = codonINFO[9] + 1;
			break;
		case	 "AUC":
			codonINFO[10] = codonINFO[10] + 1;
			break;
		case	 "AUA":
			codonINFO[11] = codonINFO[11] + 1;
			break;
		case	 "AUG":
			codonINFO[12] = codonINFO[12] + 1;
			break;
		case	 "GUU":
			codonINFO[13] = codonINFO[13] + 1;
			break;
		case	 "GUC":
			codonINFO[14] = codonINFO[14] + 1;
			break;
		case	 "GUA":
			codonINFO[15] = codonINFO[15] + 1;
			break;
		case	 "GUG":
			codonINFO[16] = codonINFO[16] + 1;
			break;
		case	 "UCU":
			codonINFO[17] = codonINFO[17] + 1;
			break;
		case	 "UCC":
			codonINFO[18] = codonINFO[18] + 1;
			break;
		case	 "UCA":
			codonINFO[19] = codonINFO[19] + 1;
			break;
		case	 "UCG":
			codonINFO[20] = codonINFO[20] + 1;
			break;
		case	 "AGU":
			codonINFO[21] = codonINFO[21] + 1;
			break;
		case	 "AGC":
			codonINFO[22] = codonINFO[22] + 1;
			break;
		case	 "CCU":
			codonINFO[23] = codonINFO[23] + 1;
			break;
		case	 "CCC":
			codonINFO[24] = codonINFO[24] + 1;
			break;
		case	 "CCA":
			codonINFO[25] = codonINFO[25] + 1;
			break;
		case	 "CCG":
			codonINFO[26] = codonINFO[26] + 1;
			break;
		case	 "ACU":
			codonINFO[27] = codonINFO[27] + 1;
			break;
		case	 "ACC":
			codonINFO[28] = codonINFO[28] + 1;
			break;
		case	 "ACA":
			codonINFO[29] = codonINFO[29] + 1;
			break;
		case	 "ACG":
			codonINFO[30] = codonINFO[30] + 1;
			break;
		case	 "UAU":
			codonINFO[31] = codonINFO[31] + 1;
			break;
		case	 "UAC":
			codonINFO[32] = codonINFO[32] + 1;
			break;
		case	 "GCU":
			codonINFO[33] = codonINFO[33] + 1;
			break;
		case	 "GCC":
			codonINFO[34] = codonINFO[34] + 1;
			break;
		case	 "GCA":
			codonINFO[35] = codonINFO[35] + 1;
			break;
		case	 "GCG":
			codonINFO[36] = codonINFO[36] + 1;
			break;
		case	 "CAU":
			codonINFO[37] = codonINFO[37] + 1;
			break;
		case	 "CAC":
			codonINFO[38] = codonINFO[38] + 1;
			break;
		case	 "CAA":
			codonINFO[39] = codonINFO[39] + 1;
			break;
		case	 "CAG":
			codonINFO[40] = codonINFO[40] + 1;
			break;
		case	 "AAU":
			codonINFO[41] = codonINFO[41] + 1;
			break;
		case	 "AAC":
			codonINFO[42] = codonINFO[42] + 1;
			break;
		case	 "AAA":
			codonINFO[43] = codonINFO[43] + 1;
			break;
		case	 "AAG":
			codonINFO[44] = codonINFO[44] + 1;
			break;
		case	 "GAU":
			codonINFO[45] = codonINFO[45] + 1;
			break;
		case	 "GAC":
			codonINFO[46] = codonINFO[46] + 1;
			break;
		case	 "GAA":
			codonINFO[47] = codonINFO[47] + 1;
			break;
		case	 "GAG":
			codonINFO[48] = codonINFO[48] + 1;
			break;
		case	 "UGU":
			codonINFO[49] = codonINFO[49] + 1;
			break;
		case	 "UGC":
			codonINFO[50] = codonINFO[50] + 1;
			break;
		case	 "UGG":
			codonINFO[51] = codonINFO[51] + 1;
			break;
		case	 "CGU":
			codonINFO[52] = codonINFO[52] + 1;
			break;
		case	 "CGC":
			codonINFO[53] = codonINFO[53] + 1;
			break;
		case	 "CGA":
			codonINFO[54] = codonINFO[54] + 1;
			break;
		case	 "CGG":
			codonINFO[55] = codonINFO[55] + 1;
			break;
		case	 "AGA":
			codonINFO[56] = codonINFO[56] + 1;
			break;
		case	 "AGG":
			codonINFO[57] = codonINFO[57] + 1;
			break;
		case	 "GGU":
			codonINFO[58] = codonINFO[58] + 1;
			break;
		case	 "GGC":
			codonINFO[59] = codonINFO[59] + 1;
			break;
		case	 "GGA":
			codonINFO[60] = codonINFO[60] + 1;
			break;
		case	 "GGG":
			codonINFO[61] = codonINFO[61] + 1;
			break;
		case	 "UAA":
			codonINFO[62] = codonINFO[62] + 1;
			break;
		case	 "UAG":
			codonINFO[63] = codonINFO[63] + 1;
			break;
		case	 "UGA":
			codonINFO[64] = codonINFO[64] + 1;
			break;
		default:
			codonINFO[0] = codonINFO[0] + 1;
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
		info = "Index" + "\t" + "Headinfo" + "\t" + "Status" + "\t" + 
				"UUU_(F)\t" + 
				"UUC_(F)\t" + 
				"UUA_(L)\t" + 
				"UUG_(L)\t" + 
				"CUU_(L)\t" + 
				"CUC_(L)\t" + 
				"CUA_(L)\t" + 
				"CUG_(L)\t" + 
				"AUU_(I)\t" + 
				"AUC_(I)\t" + 
				"AUA_(I)\t" + 
				"AUG_(M)\t" + 
				"GUU_(V)\t" + 
				"GUC_(V)\t" + 
				"GUA_(V)\t" + 
				"GUG_(V)\t" + 
				"UCU_(S)\t" + 
				"UCC_(S)\t" + 
				"UCA_(S)\t" + 
				"UCG_(S)\t" + 
				"AGU_(S)\t" + 
				"AGC_(S)\t" + 
				"CCU_(P)\t" + 
				"CCC_(P)\t" + 
				"CCA_(P)\t" + 
				"CCG_(P)\t" + 
				"ACU_(U)\t" + 
				"ACC_(U)\t" + 
				"ACA_(U)\t" + 
				"ACG_(U)\t" + 
				"UAU_(Y)\t" + 
				"UAC_(Y)\t" + 
				"GCU_(A)\t" + 
				"GCC_(A)\t" + 
				"GCA_(A)\t" + 
				"GCG_(A)\t" + 
				"CAU_(H)\t" + 
				"CAC_(H)\t" + 
				"CAA_(Q)\t" + 
				"CAG_(Q)\t" + 
				"AAU_(N)\t" + 
				"AAC_(N)\t" + 
				"AAA_(K)\t" + 
				"AAG_(K)\t" + 
				"GAU_(D)\t" + 
				"GAC_(D)\t" + 
				"GAA_(E)\t" + 
				"GAG_(E)\t" + 
				"UGU_(C)\t" + 
				"UGC_(C)\t" + 
				"UGG_(W)\t" + 
				"CGU_(R)\t" + 
				"CGC_(R)\t" + 
				"CGA_(R)\t" + 
				"CGG_(R)\t" + 
				"AGA_(R)\t" + 
				"AGG_(R)\t" + 
				"GGU_(G)\t" + 
				"GGC_(G)\t" + 
				"GGA_(G)\t" + 
				"GGG_(G)\t" + 
				"UAA_(N/A)\t" + 
				"UAG_(N/A)\t" + 
				"UGA_(N/A)\n";
		write(outPath,info);
		readmRNA(mRNAPath);
	}
}
