import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class GOEntry {
	
	private String id;
	private String name;
	private HashSet<String> isa;
	private HashSet<String> pred;
	private HashSet<String> succ;
	private HashSet<String> gene;
	
	private HashMap<String,Integer> dist;
	
	private int size;
	private int noverlap;
	private boolean istrue;
	
	private HashMap<String,ArrayList<String>> paths;
	
	public GOEntry(String id, String name, HashSet<String> isa){
		this.id = id;
		this.name = name;
		this.isa = isa;

		pred = new HashSet<String>();
		succ = new HashSet<String>();
		gene = new HashSet<String>();
		
		dist = new HashMap<String,Integer>();
		
		paths = new HashMap<String,ArrayList<String>>();
		
	}
	
	public String getId(){
		return id;
	}
	
	public String getName(){
		return name;
	}
	
	public HashSet<String> is_a(){
		return isa;
	}
	
	public HashSet<String> getPred(){
		return pred;
	}
	
	public HashSet<String> getSucc(){
		return succ;
	}
	
	public HashSet<String> getGenes(){
		return gene;
	}
	
	public void addGene(String gene){
		this.gene.add(gene);
	}
	
	public void addSucc(HashSet<String> succs){
		succ.addAll(succs);
	}
	
	public void addSucc(String succ){
		this.succ.add(succ);
	}
	
	public void addPred(HashSet<String> preds){
		pred.addAll(preds);
	}
	
	public void addPred(String pred){
		this.pred.add(pred);
	}
	
	public String toString(){
		String brk ="\n";
		String tab = "\t";
		StringBuilder entryBuilder = new StringBuilder("Id:");
		entryBuilder.append(tab);
		entryBuilder.append(id);
		entryBuilder.append(brk);
		entryBuilder.append("Name:");
		entryBuilder.append(tab);
		entryBuilder.append(name);
		entryBuilder.append(brk);
		entryBuilder.append("is a:");
		entryBuilder.append(tab);
		entryBuilder.append(isa.toString());
		entryBuilder.append(brk);
		entryBuilder.append("Pred:");
		entryBuilder.append(tab);
		entryBuilder.append(pred.toString());
		entryBuilder.append(brk);
		entryBuilder.append("Succ:");
		entryBuilder.append(tab);
		entryBuilder.append(succ.toString());
		entryBuilder.append(brk);
		entryBuilder.append("Genes:");
		entryBuilder.append(tab);
		entryBuilder.append(gene.toString());
		entryBuilder.append(brk);
		entryBuilder.append("Dist:");
		entryBuilder.append(tab);
		entryBuilder.append(dist.toString());
		entryBuilder.append(brk);
		entryBuilder.append("Size:");
		entryBuilder.append(tab);
		entryBuilder.append(this.getSize());
		
		return entryBuilder.toString();
	}
	
	public HashMap<String,Integer> getDistance(){
		return dist;
	}
	
	public void addDistance(HashMap<String,Integer> parentDist){
		int value;
		for(String pred: parentDist.keySet()){
			value = parentDist.get(pred)+1;
			if(dist.containsKey(pred)){
				if(dist.get(pred) > value){
					dist.put(pred, value);
				}
			}
			else{
				dist.put(pred, value);
			}
		}
	}
	
	public void putDistance(String id,int dist){
		this.dist.put(id,dist);
	}
	
	public void setSize(int size){
		this.size = size;
	}
	
	public int getSize(){
		return size;
	}
	
	public void setNOverlap(int noverlap){
		this.noverlap = noverlap;
	}
	
	public int getNOverlap(){
		return noverlap;
	}
	
	public void setTruth(boolean truth){
		istrue = truth;
	}
	
	public boolean getTruth(){
		return istrue;
	}
	
	public void setPath(ArrayList<String> path){
		
//		System.out.println("SetPath: " + id);
//		System.out.println(path.toString());
		for(int i =path.size()-1; i >= 0; i--){
			paths.put(path.get(i),new ArrayList<String>(path.subList(0, i+1)));
		}
//		System.out.println(paths.toString());
	}
	
	public ArrayList<String> getPath(String id){
		if(paths.containsKey(id)){
			return paths.get(id);
		}
		else{
			return null;
		}
	}
	
	public HashMap<String,ArrayList<String>> getPaths(){
		return paths;
	}

}
