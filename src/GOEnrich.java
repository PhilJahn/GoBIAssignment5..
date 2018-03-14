import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.zip.*;

public class GOEnrich {
	public enum Type{ensembl,go};
	
	public static void main(String[] args) {
		
		String oboPath ="";
		String goname ="";
		String mappingPath ="";
		String simulationPath = "";
		
		String outputPath = "";
		String ooutputPath = "";
		
		int minsize = -1;
		int maxsize = -1;
		
		Type type = null;

		
		for(int i =0; i < args.length-1; i++){
			if(args[i].equals("-obo")){
				oboPath = args[i+1];
				i++;
			}
			else if(args[i].equals("-root")){
				goname = args[i+1];
				i++; 
			}
			else if(args[i].equals("-mapping")){
				mappingPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-mappingtype")){
				if("ensembl".equals(args[i+1])){
					type = Type.ensembl;
				}
				else if("go".equals(args[i+1])){
					type = Type.go;
				}
				i++; 
			}
			else if(args[i].equals("-enrich")){
				simulationPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-o")){
				outputPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-overlapout")){
				ooutputPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-minsize")){
				minsize = Integer.parseInt(args[i+1]);
				i++; 
			}
			else if(args[i].equals("-maxsize")){
				maxsize = Integer.parseInt(args[i+1]);
				i++; 
			}
		}
		
//		System.out.println(oboPath);
//		System.out.println(goname);
//		System.out.println(mappingPath);
//		System.out.println(simulationPath);
//		System.out.println(outputPath);
//		System.out.println(minsize);
//		System.out.println(maxsize);
		
		if(oboPath.equals("") || goname.equals("") || mappingPath.equals("") || simulationPath.equals("") || outputPath.equals("") || minsize == -1 || maxsize == -1){
			System.out.println("Usage Info:\n-obo <obo_file>\n-root <GO namespace>\n-mapping <gene2go_mapping>\n-mappingtype [ensembl|go]\n[-overlapout overlap_out_tsv]\n-enrich <simulation file>\n-o <output_tsv>\n-minsize <int>\n-maxsize <int>");
		}
		else{
			
			Path oboFilePath = Paths.get(oboPath);
			
			File oboFile = oboFilePath.toFile();
			try {
				GOEnrich goe = new GOEnrich(oboFile, goname, type, mappingPath);
				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	HashMap<String,GOEntry> go_entries;
	HashMap<String,HashSet<String>> overlaps;
	HashMap<String,Gene> genes;
	
	public GOEnrich(File oboFile, String goname, Type mtype, String mappingPath) throws IOException{
		
		go_entries = new HashMap<String,GOEntry>();
		overlaps = new HashMap<String,HashSet<String>>();
		genes = new HashMap<String,Gene>();
		
		BufferedReader obobr = new BufferedReader (new FileReader(oboFile));
		String line ="";
		boolean inGOTerm = true;
		HashSet<String> goterm = new HashSet<String>();
    	String go_id = "";
    	String go_name = "";
    	String go_namespace ="";
    	HashSet<String> isa = null;
	    while ((line = obobr.readLine()) != null){
	    	inGOTerm = !(line.contains("[T"));
	        if(inGOTerm){
	        	goterm.add(line);
	        }
	        else{
	        	isa = new HashSet<String>();
	        	for(String ti: goterm){
	        		if(ti.contains("id:")){
	        			String[] tisplit = ti.split(": ");
	        			go_id = tisplit[1];
	        		}
	        		else if(ti.contains("name:")){
	        			String[] tisplit = ti.split(": ");
	        			go_name = tisplit[1];
	        		}
	        		else if(ti.contains("namespace:")){
	        			String[] tisplit = ti.split(": ");
	        			go_namespace = tisplit[1];
	        		}
	        		else if(ti.contains("is_a:")){
	        			String[] tisplit = ti.split(": ");
	        			tisplit = tisplit[1].split(" ! ");
	        			isa.add(tisplit[0]);
	        		}
	        	}
	        	if(!go_id.equals("")){
	        		if(go_namespace.equals(goname)){
	        			GOEntry newGOEntry = new GOEntry(go_id,go_name,isa);
	        			go_entries.put(go_id,newGOEntry);
	        		}
	        	}
	        	
	        	go_id = "";
	        	go_name = "";
	        	go_namespace ="";
	        	goterm.clear();
	        }
	    }
	    
	    obobr.close();
	    
	    for(String go_entry: go_entries.keySet()){
	    	setGOEntry(go_entry);
	    }
	    line ="";
	    if(mtype==Type.ensembl){
			Path mappingFilePath = Paths.get(mappingPath);
			
			File mappingFile = mappingFilePath.toFile();
			BufferedReader mapbr = new BufferedReader (new FileReader(mappingFile));
			line = mapbr.readLine();
			HashSet<String> curGOEntry = new HashSet<String>();
			while ((line = mapbr.readLine()) != null){
				String[] lineSplit = line.split("\t");
				String geneid = lineSplit[1];
				if(!geneid.equals("")){
					Gene curGene = new Gene(geneid);
					genes.put(geneid,curGene);
					String[] goes = lineSplit[2].split("\\|");
					for(int i = 0; i < goes.length; i++){
//						System.out.println(goes[i]);
						if(go_entries.containsKey(goes[i])){
    						GOEntry goe = go_entries.get(goes[i]);
    						goe.addGene(geneid);
    						for(String predgoeid : goe.getPred()){
    							go_entries.get(predgoeid).addGene(geneid);
    						}
    						curGOEntry.add(goes[i]);
	    				}
					}
					
				}
				for(String goe1: curGOEntry){
					for(String goe2: curGOEntry){
						addOverlap(goe1,goe2);
					}
				}
				curGOEntry.clear();
			}
			mapbr.close();
	    }
	    else{
	    	FileInputStream gafinput = new FileInputStream(mappingPath);
	    	GZIPInputStream gafzip = new GZIPInputStream(gafinput);
	    	InputStreamReader ireader = new InputStreamReader(gafzip);
	    	BufferedReader gafbr = new BufferedReader(ireader);

	    	String curGene = "";
	    	HashSet<String> curGOEntry = new HashSet<String>();
	    	while ((line = gafbr.readLine()) != null) {
	    		if(!line.startsWith("!")){
	    			String[] lineSplit = line.split("\t");
	    			if(lineSplit[3].equals("")){
	    				if(!curGene.equals(lineSplit[2])){
	    					Gene newGene = new Gene(curGene);
	    					genes.put(curGene,newGene);
	    					for(String goeid: curGOEntry){
	    						GOEntry goe = go_entries.get(goeid);
	    						goe.addGene(curGene);
	    						for(String predgoeid : goe.getPred()){
	    							go_entries.get(predgoeid).addGene(curGene);
	    						}
	    					}
	    					for(String goe1: curGOEntry){
	    						for(String goe2: curGOEntry){
	    							addOverlap(goe1,goe2);
	    						}
	    					}
	    					curGene = "";
	    					curGOEntry.clear();
	    				}
	    				curGene = lineSplit[2];
	    				if(go_entries.containsKey(lineSplit[4])){
	    					curGOEntry.add(lineSplit[4]);
	    				}
	    			}
	    		}
	    	}
	    	Gene newGene = new Gene(curGene);
			genes.put(curGene,newGene);
			for(String goeid: curGOEntry){
				GOEntry goe = go_entries.get(goeid);
				goe.addGene(curGene);
				for(String predgoeid : goe.getPred()){
					go_entries.get(predgoeid).addGene(curGene);;
				}
			}
			for(String goe1: curGOEntry){
				for(String goe2: curGOEntry){
					addOverlap(goe1,goe2);
				}
			}
			curGene = "";
			curGOEntry = new HashSet<String>();
			gafbr.close();
			ireader.close();
			gafzip.close();
			gafinput.close();
	    }
	    
	    
//	    for( String go1 : overlaps.keySet()){
//	    	System.out.println(go1 + ": " + overlaps.get(go1).toString());
//	    }
	    
	    for(String go_entry: go_entries.keySet()){
	    	System.out.println(go_entries.get(go_entry).toString());
	    }
		
	}
	
	@SuppressWarnings("unchecked")
	public HashSet<String> setGOEntry(String id, HashSet<String> succs){
		GOEntry curGOE = go_entries.get(id);
		HashSet<String> parents = curGOE.is_a();
		succs.add(id);
		for(String parentId: parents){
			GOEntry parentGOE = go_entries.get(parentId);
			if(parentGOE  != null){
				parentGOE.addSucc(succs);
				HashSet<String> predIds = parentGOE.getPred();
				if(predIds.size() > 0){
					for(String predId : predIds){
						if(go_entries.get(predId)  != null){
							go_entries.get(predId).addSucc(succs);
						}
					}
					curGOE.addPred(predIds);
				}
				else{
					HashSet<String> newSuccs = (HashSet<String>) parentGOE.getSucc().clone();
					curGOE.addPred(setGOEntry(parentId, newSuccs));
				}
				curGOE.addPred(parentId);
			}
		};
		return curGOE.getPred();
	}
	
	@SuppressWarnings("unchecked")
	public void setGOEntry(String id){
		GOEntry curGOE = go_entries.get(id);
		if(curGOE.getPred().size() == 0){
		HashSet<String> parents = curGOE.is_a();
		HashSet<String> succs  = new HashSet<String>();
		succs.add(id);
		for(String parentId: parents){
			GOEntry parentGOE = go_entries.get(parentId);
			if(parentGOE  != null){
				parentGOE.addSucc(succs);
				HashSet<String> predIds = parentGOE.getPred();
				if(predIds.size() > 0){
					for(String predId : predIds){
						if(go_entries.get(predId)  != null){
							go_entries.get(predId).addSucc(succs);
						}
					}
					curGOE.addPred(predIds);
				}
				else{
					HashSet<String> newSuccs = new HashSet<String>(parentGOE.getSucc());
					curGOE.addPred(setGOEntry(parentId, newSuccs));
				}
				curGOE.addPred(parentId);
			}
		}
		}
	}
	
	public void addOverlap(String goe1, String goe2){
		boolean already_included = false;
		// goe2 before goe1
		if(goe1.compareTo(goe2) < 0){
			if(overlaps.containsKey(goe2)){
				already_included = overlaps.get(goe2).contains(goe1);
			}
		}
		else{
			if(overlaps.containsKey(goe1)){
				already_included = overlaps.get(goe1).contains(goe2);
			}
		}
		if(!already_included){
			HashSet<String> set1 = new HashSet<String>(go_entries.get(goe1).getPred());
			set1.add(goe1);
			HashSet<String> set2 = new HashSet<String>(go_entries.get(goe2).getPred());
			set2.add(goe2);
			
			for(String ov1: set1){
				for(String ov2: set2){
					if(ov1.compareTo(ov2) < 0){
						if(!overlaps.containsKey(ov2)){
							overlaps.put(ov2,new HashSet<String>());
						}
						overlaps.get(ov2).add(ov1);
					}
					else{
						if(!overlaps.containsKey(ov1)){
							overlaps.put(ov1,new HashSet<String>());
						}
						overlaps.get(ov1).add(ov2);						
					}
				}
			}
		}
		
	}
	
}
