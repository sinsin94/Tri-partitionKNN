
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Arrays;
import coser.common.SimpleTool;
import weka.clusterers.EM;
import weka.clusterers.FarthestFirst;
import weka.clusterers.HierarchicalClusterer;
import weka.clusterers.SimpleKMeans;
import weka.core.Instance;
import weka.core.Instances;
import weka.filters.unsupervised.attribute.Remove;

public class RDC extends Instances {
	private static final long serialVersionUID = -82397417233430226L;
	private static String DataSet = "sonar";
	private static String fileAddress = "/Users/dengsiyu/eclipse/workspace/Coser_Triple/data/voting.arff";
	public static int numClass;
	
	protected double[][] similarityMatrix;
	protected double[][] powerMatrix;
	protected double[][] modifiedMatrix;
	protected double[] rankVector;
	protected double zuni = 0.95;
	protected double threshold = 0.01;
	protected static double percentage = 0.2;
	protected int times = 0;
	protected double[] thresholdVector;
	protected int[] neighborCountVector;
	protected int[][] neighborMatrix;
	
	public RDC (Reader paraReader) throws Exception, IOException {
		super(paraReader);
	}// Of constructor
	
	public RDC (Instances instances) throws Exception, IOException {
		super(instances);
	}// Of constructor
	
	public static double getSimilarityValue (Instance former, Instance latter) {
		int count = 0;
		for (int i = 0; i < former.numAttributes() - 1; i ++) {
			int valueFormer = (int)former.value(i);
			int valueLatter = (int)latter.value(i);
			if (valueFormer == valueLatter) {
				count ++;
			}// Of if
		}//Of for i
		
		return (count + 0.0) / (former.numAttributes() - 1);
	}// Of getSimilarityValue
	
	public void getSimilarityMatrix () {
		similarityMatrix = new double[numInstances()][numInstances()];
		for (int i = 0; i < numInstances(); i ++) {
			for (int j = i; j < numInstances(); j ++) {
				similarityMatrix[i][j] = getSimilarityValue(instance(i), instance(j));
				similarityMatrix[j][i] = similarityMatrix[i][j];
			}// Of for j
		}// Of for i
	}// Of getSimilarityMatrix
	
	public void asignMostSimilar () {
		double[] mostSimilarVector = new double[numInstances()];
		for (int i = 0; i < numInstances(); i ++) {
			double max = Double.MIN_VALUE;
			for (int j = 0; j < numInstances(); j ++) {
				if (i != j) {
					if (similarityMatrix[i][j] > max) {
						max = similarityMatrix[i][j];
					}// Of if
				}// Of if
			}// Of for j 
			mostSimilarVector[i] = max;
		}// Of for i

		boolean[][] statuMatrix = new boolean[numInstances()][numInstances()];
		int[] countVector = new int[numInstances()];
		for (int i = 0; i < numInstances(); i ++) {			
			for (int j = 0; j < numInstances(); j ++) {
				if (similarityMatrix[i][j] >= mostSimilarVector[i]) {             				
					countVector[i] ++;
					statuMatrix[i][j] = true;
				}// Of if
			}// Of for j		
		}// Of for i
		
		powerMatrix = new double[numInstances()][numInstances()];
		for (int i = 0; i < numInstances(); i ++) {
			double value = (1 + 0.0) / countVector[i];
			for (int j = 0; j < numInstances(); j ++) {
				if (statuMatrix[i][j]) {
					powerMatrix[i][j] = value;
				}// Of if 
			}// Of for j
		}// Of for i
	}// Of asignMostSimilar
	
	public void revertPowerMatrix () {
		double[][] tempMatrix = new double[numInstances()][numInstances()];
		for (int i = 0; i < numInstances(); i ++) {
			for (int j = 0; j < numInstances(); j ++) {
				tempMatrix[i][j] = powerMatrix[j][i];
			}// Of for j
		}// Of for i
		powerMatrix = tempMatrix;
	}// Of revertPowerMatrix
	
	public void getModifiedMatrix () {
		modifiedMatrix = new double[numInstances()][numInstances()];
		double[][] EETNMatrix = new double[numInstances()][numInstances()];
		for (int i = 0; i < numInstances(); i ++) {
			for (int j = 0; j < numInstances(); j ++) {
				EETNMatrix[i][j] = (1 + 0.0) / numInstances();
			}// Of for j
		}// Of for i
		
		for (int i = 0; i < numInstances(); i ++) {
			for (int j = 0; j < numInstances(); j ++) {
				modifiedMatrix[i][j] = powerMatrix[i][j] * zuni + (1 - zuni) * EETNMatrix[i][j];
			}// Of for j
		}// Of for i
	}// Of getModifiedMatrix
	
	public static double getLevelValue (double[] vectorOne, double[] vectorTwo) {
		double sum = 0;
		for (int i = 0; i < vectorOne.length; i ++) {
			sum += Math.abs(vectorOne[i] - vectorTwo[i]);
		}// Of for i
		
		return sum;
	}// Of getLevelValue
	
	public void getRankVector () {
		rankVector = new double[numInstances()];
		for (int i = 0; i < numInstances(); i ++) {
			rankVector[i] = 1;
		}// Of for i		
		times = 0;
		double currentOrderNorm = Double.MAX_VALUE;
		while (currentOrderNorm > threshold) {
			double[] tempRankVector = new double[numInstances()];
			for (int i = 0; i < numInstances(); i ++) {
				for (int j = 0; j < numInstances(); j ++) {
					tempRankVector[i] += modifiedMatrix[i][j] * rankVector[j];
				}// Of for j
			}// Of for i
				
			currentOrderNorm = 0;
			for (int i = 0; i < rankVector.length; i ++) {
				currentOrderNorm += Math.abs(tempRankVector[i] - rankVector[i]);
			}// Of for i
				
			for (int i = 0; i < rankVector.length; i ++) {
				rankVector[i] = tempRankVector[i];
			}// Of for i
			
			times ++;
		}// Of while
	}// Of getRankVector
	
	public RP[] getRepresentativeInstance () {
		RP[] reperesentativeSet = new RP[numInstances()];
		for (int i = 0; i < numInstances(); i ++) {
			reperesentativeSet[i] = new RP();
			reperesentativeSet[i].instance = instance(i);
			reperesentativeSet[i].index = i;
			reperesentativeSet[i].rankValue = rankVector[i];
		}// Of for i
		
		return reperesentativeSet;
	}// Of getRepresentativeInstance
	
	public static RP[] getSortRepresentativeInstance (RP[] representativeSet) {
		RP tempRepresentative = new RP();
		for (int i = 0; i < representativeSet.length; i ++) {
			for (int j = i + 1; j < representativeSet.length; j ++) {
				if (representativeSet[i].rankValue <representativeSet[j].rankValue) {
					tempRepresentative = representativeSet[i];
					representativeSet[i] = representativeSet[j];
					representativeSet[j] = tempRepresentative;
				}// Of if
			}// Of for j
		}// Of for i
		
		return representativeSet;
	}// Of getSortRepresentativeInstance
	
	public int[] getCriticalInstances (double percentage) {
		RP[] representativeSet = getRepresentativeInstance();
		RP[] sortRepresentativeSet = getSortRepresentativeInstance(representativeSet);
		
		int count = (int)(percentage * sortRepresentativeSet.length);
		RP[] criticalSet = new RP[count];
		RP[] commonSet = new RP[sortRepresentativeSet.length - count];
		
		int countCritical = 0;
		int countCommon = 0;
		for (int i = 0; i < sortRepresentativeSet.length; i ++) {
			if (i < count) {
				criticalSet[countCritical] = sortRepresentativeSet[i];
				countCritical ++;
			} else {
				commonSet[countCommon] = sortRepresentativeSet[i];
				countCommon ++;
			}// Of if
		}// Of for i
		
		double[][] matrix = new double[criticalSet.length][criticalSet.length];
		for (int i = 0; i < criticalSet.length; i ++) {
			Instance tempOne = criticalSet[i].instance;
			for (int j = 0; j < criticalSet.length; j ++) {
				Instance tempTwo = criticalSet[j].instance;
				matrix[i][j] = getSimilarityValue(tempOne, tempTwo);
			}// Of for j
		}// Of for i
		
		double minValue = Double.MAX_VALUE;
		int iIndex = -1;
		int jIndex = -1;
		for (int i = 0; i < criticalSet.length; i ++) {
			for (int j = 0; j < criticalSet.length; j ++) {
				if (matrix[i][j] < minValue) {
					minValue = matrix[i][j];
					iIndex = i;
					jIndex = j;
				}// Of if
			}// Of for j
		}// Of for i
		
		Instance criticalInstanceOne = criticalSet[iIndex].instance;
		Instance criticalInstanceTwo = criticalSet[jIndex].instance;
		int uncertainCount = 0;
		for (int i = 0; i < criticalSet.length; i ++) {
			double similarityToOne = getSimilarityValue(criticalSet[i].instance, criticalInstanceOne);
			double similarityToTwo = getSimilarityValue(criticalSet[i].instance, criticalInstanceTwo);
			if (similarityToOne > similarityToTwo) {
				criticalSet[i].predictedLabel = 0;
			} else {
				criticalSet[i].predictedLabel = 1;
			}// Of if
			
			if (Math.abs(similarityToOne - similarityToTwo) < 0.2) {
				criticalSet[i].predictedLabel = -1;
				uncertainCount ++;
			}// Of if
		}// Of for i
	
		System.out.println("AA " + uncertainCount);
		
		if (uncertainCount != 0) {
			for (int i = 0; i < criticalSet.length; i ++) {
				if (criticalSet[i].predictedLabel == -1) {
					Instance instanceOne = criticalSet[i].instance;
					int countOne = 0;
					int countTwo = 0;
					double sumOne = 0;
					double sumTwo = 0;
					for (int j = 0; j < criticalSet.length; j ++) {
						if (criticalSet[j].predictedLabel == 0) {
							Instance instanceTwo = criticalSet[j].instance;
							sumOne += getSimilarityValue(instanceOne, instanceTwo);
							countOne ++;
						}// Of if
						if (criticalSet[j].predictedLabel == 1) {
							Instance instanceTwo = criticalSet[j].instance;
							sumTwo += getSimilarityValue(instanceOne, instanceTwo);
							countTwo ++;
						}// Of if
					}// Of for j
					sumOne = sumOne / countOne;
					sumTwo = sumTwo / countTwo;
					if (sumOne > sumTwo) {
						criticalSet[i].predictedLabel = 0;
					} else {
						criticalSet[i].predictedLabel = 1;
					}// Of if
				}// Of if
			}// Of for i
		}// Of if
		
		for (int i = 0; i < criticalSet.length; i ++) {
			System.out.print(criticalSet[i].predictedLabel + " ");
			System.out.println((int)criticalSet[i].instance.classValue());
		}// Of for i
		
		matrix = new double[commonSet.length][criticalSet.length];
		for (int i = 0; i < commonSet.length; i ++) {
			Instance tempInstance = commonSet[i].instance;
			double similarityToOne = 0;
			double similarityToTwo = 0;
			int countOne = 0;
			int countTwo = 0;
			for (int j = 0; j < criticalSet.length; j ++) {
				if (criticalSet[j].predictedLabel == 0) {
					similarityToOne += getSimilarityValue(tempInstance, criticalSet[j].instance);
					countOne ++;
				} else if (criticalSet[j].predictedLabel == 1){
					similarityToTwo += getSimilarityValue(tempInstance, criticalSet[j].instance);
					countTwo ++;
				}// Of if
			}// Of for j
			
			if (similarityToOne / countOne > similarityToTwo / countTwo) {
				commonSet[i].predictedLabel = 0;
			} else {
				commonSet[i].predictedLabel = 1;
			}// Of if
		}// Of for i
		
		RP[] resultSet = new RP[numInstances()];
		for (int i = 0; i < count; i ++) {
			resultSet[i] = criticalSet[i];
		}// Of for i
		for (int i = 0; i < numInstances() - count; i ++) {
			resultSet[i + count] = commonSet[i];
		}// Of for i
		
		RP tempObject = null;
		for (int i = 0; i < numInstances(); i ++) {
			for (int j = i; j < numInstances(); j ++) {
				if (resultSet[i].index > resultSet[j].index) {
					tempObject = resultSet[i];
					resultSet[i] = resultSet[j];
					resultSet[j] = tempObject;
				}// Of if
			}// Of for j
		}// Of for i
		
		int[] predictedLabelVector = new int[numInstances()];
		for (int i = 0; i < numInstances(); i ++) {
			predictedLabelVector[i] = resultSet[i].predictedLabel;
		}// Of for i
		
		return predictedLabelVector;
	}// Of getCriticalInstances
	
	public void delete(boolean[] paraIndices) {
		for (int i = numInstances() - 1; i >= 0; i--) {
			if (!paraIndices[i]) {
				delete(i);
			}// Of if
		}// Of for i
	}// Of delete
	
	public RDC[] divideInTwo(double paraPercentage) throws Exception {
		boolean[] firstInclusionArray = SimpleTool.generateBooleanArrayForDivision(numInstances(), paraPercentage);
		
		RDC firstDecisionSystem = new RDC(this);
		
		firstDecisionSystem.delete(firstInclusionArray);

		boolean[] secondInclusionArray = SimpleTool.revertBooleanArray(firstInclusionArray);
		RDC secondDecisionSystem = new RDC(this);
		secondDecisionSystem.delete(secondInclusionArray);

		RDC[] subsets = new RDC[2];
		subsets[0] = firstDecisionSystem;
		subsets[1] = secondDecisionSystem;

		return subsets;
	}// Of divideInTwo
	
	public int[] getRankDividedSets (double percentage) {
		getSimilarityMatrix();
		asignMostSimilar();
		revertPowerMatrix();
		getModifiedMatrix();
		getRankVector();
		int[] predictedLabelVector = this.getCriticalInstances(percentage);
		
		return predictedLabelVector;
	}// Of getRankDividedSets
	
	public RDC[] getRandomDividedSets (double percentage) {
		RDC[] randomDividedSets = null;
		try {
			randomDividedSets = divideInTwo(percentage);
		} catch (Exception e) {
			System.out.println(e.getStackTrace());
			System.out.println("Error occured in getRandomDividedSets(double).1");
		}// Of try
		
		return randomDividedSets;
	}// Of getRandomDividedSets
	
	public int[] getActualLabelVector () {
		
		int[] labelVector = new int[numInstances()];
		for (int i = 0; i < numInstances(); i ++) {
			labelVector[i] = (int)instance(i).classValue();
		}// Of for i
		
		return labelVector;
	}// Of getActualLabelVector
	
	public static Instances generateUnsupervisedSet (RDC dataSet) {
		Instances initialSet = new Instances(dataSet);
		Instances clusterSet = new Instances(dataSet);
		clusterSet.delete();
		try {
			Remove filter = new Remove();
            filter.setAttributeIndices("" + (initialSet.classIndex() + 1));
            filter.setInputFormat(initialSet);
            clusterSet = weka.filters.Filter.useFilter(initialSet, filter);
		} catch (Exception e) {
			System.out.println(e.getStackTrace());
			System.out.println("Error occured in generateUnsupervisedSet().1");
		}// Of try
	
		return clusterSet;
	}// Of generateUnsupervisedSet
	
	public static void experiment () {
		RDC dataSet = null;
		try {
			FileReader fileReader = new FileReader(fileAddress);
			dataSet = new RDC(fileReader);
			fileReader.close();
			
			dataSet.setClassIndex(dataSet.numAttributes() - 1);
		} catch (Exception e) {
			System.out.println("Error occured in experiment().1!");
			e.getStackTrace();
		}// Of try
		
		RDC[] tempSets = null;
		RDC usedSet = null;
		try {
			tempSets = dataSet.divideInTwo(0.99);
			usedSet = tempSets[0];
		} catch (Exception e) {
			System.out.println("Error occured in experiment().2!");
			e.getStackTrace();
		}// OF try
		
		
		numClass = usedSet.attribute(dataSet.numAttributes() - 1).numValues();
		
		int[] actualLabelVector = usedSet.getActualLabelVector();
		int[] rankPredictedLabelVector = usedSet.getRankDividedSets(percentage);
		Instances clusterSet = generateUnsupervisedSet(usedSet);
		int[] KMLabelVector = simpleKMeansCluster(clusterSet);
		int[] EMLabelVector = EMCluster(clusterSet);
		int[] FFLabelVector = farthestFirstCluster(clusterSet);
		int[] HCLabelVector = hierachicalCluster(clusterSet);
		
		double[][] evaluationVector = new double[5][];
		evaluationVector[0] = getEvaluation(actualLabelVector, rankPredictedLabelVector);
		evaluationVector[1] = getEvaluation(actualLabelVector, KMLabelVector);
		evaluationVector[2] = getEvaluation(actualLabelVector, EMLabelVector);
		evaluationVector[3] = getEvaluation(actualLabelVector, FFLabelVector);
		evaluationVector[4] = getEvaluation(actualLabelVector, HCLabelVector);
		
		for (int i = 0; i < 5; i ++) {
			System.out.println(Arrays.toString(evaluationVector[i]));	
		}// Of for i
		
	}// Of experiment

	public static double[] getEvaluation (int[] actualVector, int[] predictedVector) {
		double[] evaluationVector = new double[3];
		int SS = 0;
		int SD = 0;
		int DS = 0;
		int DD = 0;
		for (int i = 0; i < actualVector.length; i ++) {
			for (int j = i + 1; j < actualVector.length; j ++) {
				if (actualVector[i] == actualVector[j]) {
					if (predictedVector[i] == predictedVector[j]) {
						SS ++;
					} else {
						SD ++;
					}// Of if
				} else {
					if (predictedVector[i] == predictedVector[j]) {
						DS ++;
					} else {
						DD ++;
					}// Of if
				}// Of if
			}// Of for j
		}// Of for i
		
		int length = actualVector.length;
		double JA = (SS + 0.0) / (SS + SD + DS);
		double FM = Math.sqrt(((SS + 0.0) / (SS + SD)) * ((SS + 0.0) / (SS + DS)));
		double RA = 2 * (SS + DD + 0.0) / (length * (length - 1));
		
		evaluationVector[0] = JA;
		evaluationVector[1] = FM;
		evaluationVector[2] = RA;
		
		return evaluationVector;
	}// Of getEvaluation
	
	public static void main (String[] args) {
		
		experiment();
		
	}// Of main
	
	public static int[] simpleKMeansCluster (Instances clusterSet) {
		int[] clusterLabel = new int[clusterSet.numInstances()];
		
		try {
			SimpleKMeans cluster = new SimpleKMeans();
			cluster.setNumClusters(numClass);
			cluster.buildClusterer(clusterSet);
			
			for (int i = 0; i < clusterSet.numInstances(); i ++) {
				clusterLabel[i] = cluster.clusterInstance(clusterSet.instance(i));
			}// Of for i
		} catch (Exception e) {
			System.out.println("Something Error Occured in simpleKMeansCluster(Instances)!");
			e.getStackTrace();
		}// Of try
		
		return clusterLabel;
	}// Of simpleKMeansCluster
	
	public static int[] EMCluster (Instances clusterSet) {
		int[] clusterLabel = new int[clusterSet.numInstances()];
		
		try {
            String[] options = weka.core.Utils.splitOptions("-I 100 -N "+ numClass +" -M 1.0E-6 -S 100");
			EM cluster = new EM();
			cluster.setOptions(options);
			cluster.buildClusterer(clusterSet);
			
			for (int i = 0; i < clusterSet.numInstances(); i ++) {
				clusterLabel[i] = cluster.clusterInstance(clusterSet.instance(i));
			}// Of for i		
		} catch (Exception e) {
			System.out.println("Something Error Occured in EMCluster(Instaces)!");
			e.getStackTrace();
		}// Of try
		
		return clusterLabel;
	}// Of EMCluster
	
	public static int[] farthestFirstCluster (Instances clusterSet) {
		int[] clusterLabel = new int[clusterSet.numInstances()];
		
		try {
            String[] options = new String[2];
            options[0] = "-S";
            options[1] = "100";
            FarthestFirst cluster = new FarthestFirst();
            cluster.setOptions(options);
            cluster.setNumClusters(numClass);
            cluster.buildClusterer(clusterSet);
            
			for (int i = 0; i < clusterSet.numInstances(); i ++) {
				clusterLabel[i] = cluster.clusterInstance(clusterSet.instance(i));
			}// Of for i
		} catch (Exception e) {
			System.out.println("Something Error Occured in farthestFirstCluster(Instances)!");
			e.getStackTrace();
		}// Of try
		
		return clusterLabel;
	}// Of farthestFirstCluster
	
	public static int[] hierachicalCluster (Instances clusterSet) {
		int[] clusterLabel = new int[clusterSet.numInstances()];
		
		try {
            String[] options = new String[2];
            options[0] = "-L";
            options[1] = "WARD";
            HierarchicalClusterer cluster = new HierarchicalClusterer();
            cluster.setOptions(options);
            cluster.setNumClusters(numClass);
            cluster.buildClusterer(clusterSet);
					
			for (int i = 0; i < clusterSet.numInstances(); i ++) {
				clusterLabel[i] = cluster.clusterInstance(clusterSet.instance(i));	
			}// Of for i
		} catch (Exception e) {
			System.out.println("Something Error Occured in hierahicalCluster(Instances)!");
			e.getStackTrace();
		}// Of try
		
		return clusterLabel;
	}// Of hierachicalCluster
}// Of Class RDC

class RP {
	public Instance instance;
	public int index;
	public double rankValue;
	public int predictedLabel;
}// Of Class RP