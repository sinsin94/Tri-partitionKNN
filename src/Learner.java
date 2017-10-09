/*
package coser.datamodel.decisionsystem;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Reader;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Random;
import java.util.function.DoublePredicate;
import java.util.function.IntPredicate;

import javax.xml.transform.Templates;

import coser.common.SimpleTool;
import weka.classifiers.CostMatrix;
import weka.classifiers.meta.MetaCost;
import weka.classifiers.trees.Id3;
import weka.classifiers.trees.J48;
import weka.core.*;

public class Learner extends CostSensitiveDecisionSystem {
	private static final long serialVersionUID = -82097437233430226L;
	public static double percentage;
	public static double stepLength;
	
	public static int noneTotalTaught;
	
	public double teachCost;
	public double misCost;
	public int dataSetType;
	public int kValue;
	public int boughtValue;
	
	private int currentBestLabel;
	private int totalGuessed;
	private int totalTaught;
	private int currentGuessed;
	private int currentTaught;
	
	public double[][] misCostMatrix;
	public int[] instanceTypes;
	public double[] expectedCosts;
	public double[][] distanceMatrix;
	
	public CostMatrix costMatrixWeka;
		
	public static final int TAUGHT = 0;
	public static final int GUESSED = 1;
	public static final int DELAYED = 2;
	
	
	public Learner(Instances paraInstances) throws IOException, Exception {
		super(paraInstances);
	}// Of the first constructor

	public Learner(Reader paraReader) throws IOException, Exception {
		super(paraReader);
	}// Of the second constructor

	public void computeBoolDistance() {
		distanceMatrix = new double[numInstances()][numInstances()];
		for (int i = 0; i < numInstances(); i ++) {
			for (int j = 0; j < numInstances(); j ++) {
				double distance = 0;
				for (int k = 0; k < numAttributes() - 1; k ++) {
					double valueA = instance(i).value(k);
					double valueB = instance(j).value(k);
					if ((int)valueA != (int)valueB) {
						distance ++;
					}// Of if
				}// Of for k
				distanceMatrix[i][j] = distance;
			}// Of for j
		}// Of for i
	}// Of computeBoolDistance
	
	public void computeEuclideanDistance() {
		distanceMatrix = new double[numInstances()][numInstances()];
		for (int i = 0; i < numInstances(); i ++) {
			for (int j = 0; j < numInstances(); j ++) {
				double distanceSquare = 0;
				for (int k = 0; k < numAttributes() - 1; k ++) {
					double valueA = instance(i).value(k);
					double valueB = instance(j).value(k);
					distanceSquare += Math.pow((valueA - valueB), 2);
				}// Of for k
				distanceMatrix[i][j] = Math.sqrt(distanceSquare);
			}// Of for j
		}// Of for i
	}// Of computeEuclideanDistance
	
	public void computeManhattanDistance() {
		distanceMatrix = new double[numInstances()][numInstances()];
		for (int i = 0; i < numInstances(); i ++) {
			for (int j = 0; j < numInstances(); j ++) {
				double distance = 0;
				for (int k = 0; k < numAttributes() - 1; k ++) {
					double valueA = instance(i).value(k);
					double valueB = instance(j).value(k);
					distance += Math.abs(valueA - valueB);
				}// Of for k
				distanceMatrix[i][j] = distance;
			}// Of for j
		}// Of for i
	}// Of computeManhattanDistance
	
	public void computeCosineDistance() {
		distanceMatrix = new double[numInstances()][numInstances()];
		for (int i = 0; i < numInstances(); i ++) {
			for (int j = 0; j < numInstances(); j ++) {
				double multiResult = 0;
				double XSquare = 0;
				double YSquare = 0;
				for (int k = 0; k < numAttributes() - 1; k ++) {
					double valueA = instance(i).value(k);
					double valueB = instance(j).value(k);
					multiResult += valueA * valueB;
					XSquare += Math.pow(valueA, 2);
					YSquare += Math.pow(valueB, 2);
				}// Of for k
				distanceMatrix[i][j] = multiResult / (Math.sqrt(XSquare) * Math.sqrt(YSquare));
			}// Of for j
		}// Of for i
	}// Of computeCosineDistance
	
	public void getDistanceMatrix() {
		switch (dataSetType) {
		case 1:
			computeBoolDistance();
			break;
		case 2:
			computeEuclideanDistance();
			break;
		case 3:
			computeManhattanDistance();
			break;
		case 4:
			computeCosineDistance();
			break;
		default :
			System.out.println("Error occured in getDistanceMatrix!");
			break;
		}// Of switch
	}// Of getDistanceMatrix
	
	public void getMisCostMatrix() {
		misCostMatrix = new double[numClasses()][numClasses()];
		for (int i = 0; i < numClasses(); i ++) {
			for (int j = 0; j < numClasses(); j ++) {
				if (i == j) {
					misCostMatrix[i][j] = 0;
				} else {
					misCostMatrix[i][j] = misCost;
				}// Of if
			}// Of for j
		}// Of for i
	}// Of getMisCostMatrix
	
	public void getInsertAttribute() {
		int classIndex = numAttributes() - 1;
		Enumeration classValuesEnumeration = attribute(classIndex).enumerateValues();
		FastVector values = new FastVector();
		while(classValuesEnumeration.hasMoreElements()) {
			values.addElement("" + classValuesEnumeration.nextElement());
		}// Of while
		String title = "Predict";
		Attribute insertAttribute = new Attribute(title, values);
		insertAttributeAt(insertAttribute, 0);
	}// Of getInsertAttribute
	
	public void getAnnounceLength() {
		instanceTypes = new int[numInstances()];
		for (int i = 0; i < numInstances(); i ++) {
			instanceTypes[i] = DELAYED; 
		}//Of for i
		
		expectedCosts = new double[numInstances()];
	}// Of getAnnounceLength
	
	public void getInitialTrainingSet(double paraPercentage) {
		boolean[] isInitialBought = null;
		try {
			isInitialBought = SimpleTool.generateBooleanArrayForDivision(numInstances(), paraPercentage);
		} catch (Exception e) {
			System.out.println("Error occured in initializeTaughtInstance!");
		}// Of try
		
		for (int i = 0; i < isInitialBought.length; i ++) {
			if (isInitialBought[i]) {
				instance(i).setValue(0, instance(i).value(numAttributes() - 1));
				instanceTypes[i] = TAUGHT;
			}// Of if
		}// Of for i
	}// Of getInitialTrainingSet
	
	public void getInitialization(double paraPercentage) {
		getDistanceMatrix();
		getMisCostMatrix();
		getInsertAttribute();
		getAnnounceLength();
		getInitialTrainingSet(paraPercentage);
	}// Of initialization
	
	public void getParameterSetting() {
		teachCost = 10;
		misCost = 50;
		dataSetType = 2;
		kValue = 15;
		stepLength = 0.02;
		
		boughtValue = (int)(stepLength * numInstances());
		
		totalTaught = 0;
		totalGuessed = 0;
	}// Of getParameterSetting
	
	public int[] getNearestNeighbors(int paraIndex) {
		int[] tempNeighbors = new int[kValue + 1];
		double[] tempDistances = new double[kValue + 1];
		int tempDeterminedNeighors = 0;
		for (int i = 0; i < numInstances(); i ++) {
			if (instanceTypes[i] == DELAYED) {
				continue;
			}//Of if
			
			double tempCurrentDistance = distanceMatrix[paraIndex][i];
			int j;
			
			for (j = tempDeterminedNeighors - 1; j >= 0; j --) {
				if (tempDistances[j] > tempCurrentDistance) {
					tempDistances[j + 1] = tempDistances[j];
					tempNeighbors[j + 1] = tempNeighbors[j];
				} else {
					break;
				}//Of if
			}//Of for j
			
			//Insert now
			tempDistances[j + 1] = tempCurrentDistance;
			tempNeighbors[j + 1] = i;
			
			if (tempDeterminedNeighors < kValue) {
				tempDeterminedNeighors ++;
			}//Of if
		}//Of for i
		
		int[] resultNeighbors =  new int[kValue];
		for (int i = 0; i < kValue; i ++) {
			resultNeighbors[i] = tempNeighbors[i];
		}// Of for i
		
		return resultNeighbors;
	}//Of getNearestNeighbors
	
	public double getExpectedCost(int[] paraNeighbors) {
		double tempCost = 0;
		double tempMinimalCost = 10000;
		int tempClassIndex = numAttributes() - 1;
		int tempClassLabels = attribute(tempClassIndex).numValues();
		int tempNeighborClass;
		
		for (int i = 0; i < tempClassLabels; i ++) {
			tempCost = 0;
			for (int j = 0; j < paraNeighbors.length; j ++) {
				tempNeighborClass = (int)instance(paraNeighbors[j]).value(tempClassIndex);
				tempCost += misCostMatrix[tempNeighborClass][i];
			}//Of for j
			
			if (tempMinimalCost > tempCost) {
				tempMinimalCost = tempCost;
				currentBestLabel = i;
			}//Of if
		}//Of for i
		
		return tempMinimalCost / kValue;
	}//Of getExpectedCost

	public void getOneRoundGuessedInstances() {
		int[] tempNeighborIndices;
		double tempMisCost;
		currentGuessed = 0;

		for (int i = 0; i < numInstances(); i ++) {
			if (instanceTypes[i] == TAUGHT) {
				continue;
			} else if (instanceTypes[i] == GUESSED){
				continue;
			}// Of if
			
			tempNeighborIndices = getNearestNeighbors(i);
			tempMisCost = getExpectedCost(tempNeighborIndices);
			expectedCosts[i] = tempMisCost;
			if (tempMisCost < teachCost) {
				instance(i).setValue(0, currentBestLabel);
				instanceTypes[i] = GUESSED;
				currentGuessed ++;
				totalGuessed ++;
			}//Of if
		}//Of for I
	}//Of getOneRoundGuessedInstances
	
	public int getOneRoundTaughtInstances() {
		int[] tempIndices = new int[boughtValue + 1];
		double[] tempCosts = new double[boughtValue + 1];
		currentTaught = 0;
		
		for (int i = 0; i < numInstances(); i ++) {
			if (instanceTypes[i] != DELAYED) {
				continue;
			}//Of if
		
			double tempExceptedMisclassificationCosts = expectedCosts[i];
			int j;
			
			for (j = currentTaught - 1; j >= 0; j--) {
				if (tempCosts[j] < tempExceptedMisclassificationCosts) {
					tempCosts[j + 1] = tempCosts[j];
					tempIndices[j + 1] = tempIndices[j];
				} else {
					break;
				}// Of if 
			}// Of for j
			
			tempCosts[j + 1] = tempExceptedMisclassificationCosts;
			tempIndices[j + 1] = i;
			
			if (currentTaught < boughtValue) {
				currentTaught ++;
			}// Of if 
		}//Of for i
		
		for (int i = 0; i < currentTaught; i ++) {
			instanceTypes[tempIndices[i]] = TAUGHT;
			//Copy the actual class label.
			instance(tempIndices[i]).setValue(0, instance(tempIndices[i]).value(numAttributes() - 1));
			totalTaught ++;
		}//Of for i
		
		return currentTaught;
	}//Of getOneRoundTaughtInstances
	
	public int[][] getIterativeLearning() {
		int[][] tempRecordMatrix = new int[numInstances() / boughtValue][2];
		int tempTimes = 0;
		
		boolean isTerminal = getTerminationSign();
		while (!isTerminal) {
			getOneRoundGuessedInstances();
			tempRecordMatrix[tempTimes][0] = currentGuessed;
			getOneRoundTaughtInstances();
			tempRecordMatrix[tempTimes][1] = currentTaught;
			isTerminal = getTerminationSign();
			tempTimes ++;
		}// Of while
		
		int[][] recordMatrix = new int[tempTimes][2];
		for (int i = 0; i < tempTimes; i ++) {
			recordMatrix[i] = tempRecordMatrix[i];
		}// Of for i
		
		return recordMatrix;
	}// Of getIterativeLearning
	
	public boolean getTerminationSign() {
		int initialInstance = (int)(percentage * numInstances());
		int total = initialInstance + totalGuessed + totalTaught;
		
		boolean flag = true;
		if(total != numInstances()) {
			flag = false;
		}// Of 
		
		return flag;
	}// Of getTerminationSign
	
	public int[] getVerifiedResult() {
		int tempCorrect = 0;
		int tempUncorrect = 0;
		int[] tempClassificationResult = new int[2];
		for (int i = 0; i < numInstances(); i ++) {
			if (instanceTypes[i] == GUESSED) {
				int tempPredictValue = (int)instance(i).value(0);
				int tempTrueValue = (int)instance(i).value(numAttributes() - 1);
				if (tempPredictValue == tempTrueValue) {
					tempCorrect ++;
				} else {
					tempUncorrect ++;
				}// Of if
			}// Of if
		}// Of for i
		
		tempClassificationResult[0] = tempCorrect;
		tempClassificationResult[1] = tempUncorrect;
		
		return tempClassificationResult;
	}// Of getVerifiedResult
	
	public static int[][][] learnerTest() {
		Learner tempLearner = null;
		Learner firstSet = null;
		Learner secondSet = null;
		Learner thirdSet = null;
		
		String arffFilename = "E:/data/wine.arff";
		
		try {
			FileReader fileReader = new FileReader(arffFilename);
			tempLearner = new Learner(fileReader);
			fileReader.close();
			
			firstSet = new Learner(tempLearner);
			secondSet = new Learner(tempLearner);
			thirdSet = new Learner(tempLearner);
		} catch (Exception ee) {
			System.out.println("Error occured in learnerTest!");
		}// Of try
		
		int[][] noneResult = compareNone(firstSet);
		int[][] kNNResult = compareKNN(secondSet);	
		int[][] TreeResult = compareDecisionTree(thirdSet);
		
		int[][][] printResult = new int[3][][];
		printResult[0] = kNNResult;
		printResult[1] = TreeResult;
		printResult[2] = noneResult;
		
		return printResult;
	}// Of labelTest

	public static void main(String[] args) {
		try {
			FileOutputStream fs = new FileOutputStream(new File("D:\\wine15.txt"));
			PrintStream p = new PrintStream(fs);
			
			for (int i = 0; i < 10; i ++) {
				percentage = (i + 1) * 0.01;
				p.println(percentage);
				for (int j = 0; j < 10; j ++) {
					int[][][] tempCubic = learnerTest();
					p.println(Arrays.deepToString(tempCubic));
				}// Of for j
			}// Of for i
			
			p.close();
		} catch (Exception e) {
			System.out.println("Error occured in main!");
			e.printStackTrace();
		}// Of try
		
		System.out.println("Finish!");
	}// Of main

	public static int[][] compareNone(Learner tempLearner) {
		tempLearner.getParameterSetting();
		tempLearner.getInitialization(percentage);
		int[][] recordMatrix = tempLearner.getIterativeLearning();
		int[] correctAndUncorrect = tempLearner.getVerifiedResult();
		
		int[][] returnResult = new int[recordMatrix.length + 2][2];
		returnResult[0] = correctAndUncorrect;
		returnResult[1][0] = tempLearner.totalGuessed;
		returnResult[1][1] = tempLearner.totalTaught;
		for (int i = 2; i < returnResult.length; i ++) {
			returnResult[i] = recordMatrix[i - 2];
		}// Of for i
		
		noneTotalTaught = tempLearner.totalTaught;
		
		return returnResult;
	}// Of compareSelf
	
	public int getNoneClassLabel(int paraIndex) {
		double[] distanceVector = distanceMatrix[paraIndex];
		
		int count = 0;
		for (int i = 0; i < numInstances(); i ++) {
			if (instanceTypes[i] == TAUGHT) {
				count ++;
			}// Of if
		}// Of for i
		
		int index = 0;
		double[][] usefulInformation = new double[count][2];
		for (int i = 0; i < numInstances(); i ++) {
			if (instanceTypes[i] == TAUGHT) {
				usefulInformation[index][0] = i;
				usefulInformation[index][1] = distanceVector[i];
				index ++;
			}// Of if
		}// Of for i
		
		for (int i = 0; i < usefulInformation.length; i ++) {
			for (int j = i + 1; j < usefulInformation.length; j ++) {
				if (usefulInformation[i][1] > usefulInformation[j][1]) {
					double[] tempArray = usefulInformation[i];
					usefulInformation[i] = usefulInformation[j];
					usefulInformation[j] = tempArray;
				}// Of if
			}// Of for j
		}// Of for i
		
		int[] neighborIndices = new int[kValue];
		int[] vote = new int[numClasses()];
		for (int i = 0; i < neighborIndices.length; i ++) {
			int classLabel = (int)instance(neighborIndices[i]).classValue();
			vote[classLabel] ++;
		}// Of for i
		
		int max = Integer.MIN_VALUE;
		int predictedLabel = -1;
		for (int i = 0; i < vote.length; i ++) {
			if (vote[i] > max) {
				max = vote[i];
				predictedLabel = i;
			}// Of if
		}// Of for i
		
		return predictedLabel;
	}// Of getNoneClassLabel
	
	public void getKNNClassification() {
		for (int i = 0; i < numInstances(); i ++) {
			if (instanceTypes[i] == TAUGHT) {
				continue;
			}// Of if
			
			int predictedLabel = getNoneClassLabel(i);
			instance(i).setValue(0, predictedLabel);
		}// Of for i
	}// Of getKNNClassification
	
	public int[] getNoneVerifiedResult() {
		int tempCorrect = 0;
		int tempUncorrect = 0;
		
		for (int i = 0; i < numInstances(); i ++) {
			int predictedLabel = (int)instance(i).value(0);
			int intrinsicLabel = (int)instance(i).classValue();
			if (predictedLabel == intrinsicLabel) {
				tempCorrect ++;
			} else if (predictedLabel != intrinsicLabel) {
				tempUncorrect ++;
			} else {
				System.out.println("Error occured in getNoneVerifiedResult!");
			}// Of if
		}// Of for i

		int[] tempClassificationResult = new int[2];
		tempClassificationResult[0] = tempCorrect;
		tempClassificationResult[1] = tempUncorrect;
		
		return tempClassificationResult;
	}// Of getNoneVerifiedResult
	
	public static int[][] compareKNN(Learner tempLearner) {
		int labeled = (int)(percentage * tempLearner.numInstances()) + noneTotalTaught;
		double newPercentage = (labeled + 0.0) / tempLearner.numInstances();
		
		tempLearner.getParameterSetting();
		tempLearner.getInitialization(newPercentage);
		tempLearner.getKNNClassification();
		int[] classificationResult = tempLearner.getNoneVerifiedResult();
		int[][] returnResult = new int[1][];
		returnResult[0] = classificationResult;
		
		return returnResult;
	}// Of compareKNN
	
	public Learner[] getDividedDataSets(double paraPercentage) throws Exception {
		boolean[] firstInclusionArray = SimpleTool.generateBooleanArrayForDivision(numInstances(), paraPercentage);
		Learner firstDecisionSystem = new Learner(this);
		firstDecisionSystem.delete(firstInclusionArray);

		boolean[] secondInclusionArray = SimpleTool.revertBooleanArray(firstInclusionArray);
		Learner secondDecisionSystem = new Learner(this);
		secondDecisionSystem.delete(secondInclusionArray);

		Learner[] subsets = new Learner[2];
		subsets[0] = firstDecisionSystem;
		subsets[1] = secondDecisionSystem;
		
		return subsets;
	}// Of getDividedDataSets
	
	public static int[][] compareDecisionTree(Learner tempLearner) {
		int labeled = (int)(percentage * tempLearner.numInstances()) + noneTotalTaught;
		double newPercentage = (labeled + 0.0) / tempLearner.numInstances();
		
		Learner trainingSet = null;
		Learner testingSet = null;
		Learner[] tempSets = null;
		
		Learner trainingCSId3 = null;
		Learner trainingCSJ48 = null;
		Learner testingCSId3 = null;
		Learner testingCSJ48 = null;
		
		try {
			tempSets = tempLearner.getDividedDataSets(newPercentage);
			trainingSet = tempSets[0];
			testingSet = tempSets[1];
			
			trainingCSId3 = new Learner(trainingSet);
			trainingCSJ48 = new Learner(trainingSet);
			testingCSId3 = new Learner(testingSet);
			testingCSJ48 = new Learner(testingSet);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error occured in compareDecisionTree!");
		}// Of try
		
		int[] CSId3Result = compareCSId3(trainingCSId3, testingCSId3);
		int[] CSJ48Result = compareCSJ48(trainingCSJ48, testingCSJ48);
		//int[][] returnResult = new int[2][];
		//returnResult[0] = CSId3Result;
		//returnResult[1] = CSJ48Result;
		int[][] returnResult = new int[2][];
		returnResult[0] = CSId3Result;
		returnResult[1] = CSJ48Result;
		
		return returnResult;
	}// Of compareDecisionTrees
	
	public Learner getCSId3TrainingSet() {
		Learner CSTrainingSet = null;
		try {
			CSTrainingSet = new Learner(this);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error occured in getCSId3TrainingSet!");
		}// Of try
		
		CSTrainingSet.deleteWithMissingClass();
		Id3 classifier = new Id3();
		
		costMatrixWeka = new CostMatrix(CSTrainingSet.numClasses());
		for (int i = 0; i < numClasses(); i ++) {
			for (int j = 0; j < numClasses(); j ++) {
				if (i == j) {
					costMatrixWeka.setCell(i, j, (0 + 0.0));
				} else {
					costMatrixWeka.setCell(i, j, (50 + 0.0));
				}// Of if
			}// Of for j
		}// Of for i
	
		Random random = null;
		if (!(classifier instanceof WeightedInstancesHandler)) {
			random = new Random(1);
		}// Of if
		try {
			Instances CSDataSet = costMatrixWeka.applyCostMatrix(CSTrainingSet, random);
			CSTrainingSet = new Learner(CSDataSet);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error occured in getCSId3TrainingSet!");
		}// Of try
		
		return CSTrainingSet;
	}// Of getCSId3TrainingSet
	
	public static int[] compareCSId3(Learner trainingSet, Learner testingSet) {
		Learner CSTrainingSet = trainingSet.getCSId3TrainingSet();
		
		int tempCorrect = 0;
		int tempUncorrect = 0;
		int tempUnclassified = 0;
		
		Id3 classifier = new Id3();
		try {
			classifier.buildClassifier(CSTrainingSet);
			for (int i = 0; i < testingSet.numInstances(); i ++) {
				int predictedLabel = (int)classifier.classifyInstance(testingSet.instance(i));
				int trueLabel = (int)testingSet.instance(i).classValue();
				if (trueLabel < 0 || trueLabel >= testingSet.numClasses()) {
					tempUnclassified ++;
				} else {
					if (trueLabel == predictedLabel) {
						tempCorrect ++;
					} else {
						tempUncorrect ++;
					}// Of if
				}// Of if
			}// Of for i
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error occured in compareCSId3!");
		}// Of try
		
		int[] returnResult = new int[3];
		returnResult[0] = tempCorrect;
		returnResult[1] = tempUncorrect;
		returnResult[2] = tempUnclassified;
		
		return returnResult;
	}// Of compareCSId3
	
	public Learner getCSJ48TrainingSet() {
		Learner CSTrainingSet = null;
		try {
			CSTrainingSet = new Learner(this);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error occured in getCSJ48TrainingSet!");
		}// Of try
		
		CSTrainingSet.deleteWithMissingClass();
		J48 classifier = new J48();
		
		costMatrixWeka = new CostMatrix(CSTrainingSet.numClasses());
		for (int i = 0; i < numClasses(); i ++) {
			for (int j = 0; j < numClasses(); j ++) {
				if (i == j) {
					costMatrixWeka.setCell(i, j, 0);
				} else {
					costMatrixWeka.setCell(i, j, 50);
				}// Of if
			}// Of for j
		}// Of for i
	
		Random random = null;
		if (!(classifier instanceof WeightedInstancesHandler)) {
			random = new Random(1);
		}// Of if
		try {
			Instances CSDataSet = costMatrixWeka.applyCostMatrix(CSTrainingSet, random);
			CSTrainingSet = new Learner(CSDataSet);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error occured in getCSId3TrainingSet!");
		}// Of try
		
		return CSTrainingSet;
	}// Of getCSJ48TrainingSet
	
	public static int[] compareCSJ48(Learner trainingSet, Learner testingSet) {
		int tempCorrect = 0;
		int tempUncorrect = 0;
		int tempUnclassified = 0;
		
		J48 classifier = new J48();
		try {
			classifier.buildClassifier(trainingSet);
			for (int i = 0; i < testingSet.numInstances(); i ++) {
				int predictedLabel = (int)classifier.classifyInstance(testingSet.instance(i));
				int trueLabel = (int)testingSet.instance(i).classValue();
				if (trueLabel < 0 || trueLabel >= testingSet.numClasses()) {
					tempUnclassified ++;
				} else {
					if (trueLabel == predictedLabel) {
						tempCorrect ++;
					} else {
						tempUncorrect ++;
					}// Of if
				}// Of if
			}// Of for i
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error occured in compareCSJ48!");
		}// Of try
		
		int[] returnResult = new int[3];
		returnResult[0] = tempCorrect;
		returnResult[1] = tempUncorrect;
		returnResult[2] = tempUnclassified;
		
		return returnResult;
	}// Of compareCSJ48
}// Of class Learner
*/

//package coser.datamodel.decisionsystem;
//
//import java.io.File;
//import java.io.FileOutputStream;
//import java.io.FileReader;
//import java.io.IOException;
//import java.io.PrintStream;
//import java.io.Reader;
//import java.util.Arrays;
//import java.util.Enumeration;
//import java.util.Random;
//
//import coser.common.SimpleTool;
//import weka.classifiers.CostMatrix;
//import weka.classifiers.trees.Id3;
//import weka.classifiers.trees.J48;
//import weka.core.*;
//
//public class Learner extends CostSensitiveDecisionSystem {
//	private static final long serialVersionUID = -82097437233430226L;
//	public static double percentage;
//	public static double stepLength;
//	
//	public static int noneTotalTaught;
//	
//	public double teachCost;
//	public double misA;
//	public double misB;
//	public int dataSetType;
//	public int kValue;
//	public int boughtValue;
//	
//	private int currentBestLabel;
//	private int totalGuessed;
//	private int totalTaught;
//	private int currentGuessed;
//	private int currentTaught;
//	
//	public double[][] misCostMatrix;
//	public int[] instanceTypes;
//	public double[] expectedCosts;
//	public double[][] distanceMatrix;
//	
//	public CostMatrix costMatrixWeka;
//		
//	public static final int TAUGHT = 0;
//	public static final int GUESSED = 1;
//	public static final int DELAYED = 2;
//	
//	
//	public Learner(Instances paraInstances) throws IOException, Exception {
//		super(paraInstances);
//	}// Of the first constructor
//
//	public Learner(Reader paraReader) throws IOException, Exception {
//		super(paraReader);
//	}// Of the second constructor
//
//	public void computeBoolDistance() {
//		distanceMatrix = new double[numInstances()][numInstances()];
//		for (int i = 0; i < numInstances(); i ++) {
//			for (int j = 0; j < numInstances(); j ++) {
//				double distance = 0;
//				for (int k = 0; k < numAttributes() - 1; k ++) {
//					double valueA = instance(i).value(k);
//					double valueB = instance(j).value(k);
//					if ((int)valueA != (int)valueB) {
//						distance ++;
//					}// Of if
//				}// Of for k
//				distanceMatrix[i][j] = distance;
//			}// Of for j
//		}// Of for i
//	}// Of computeBoolDistance
//	
//	public void computeEuclideanDistance() {
//		distanceMatrix = new double[numInstances()][numInstances()];
//		for (int i = 0; i < numInstances(); i ++) {
//			for (int j = 0; j < numInstances(); j ++) {
//				double distanceSquare = 0;
//				for (int k = 0; k < numAttributes() - 1; k ++) {
//					double valueA = instance(i).value(k);
//					double valueB = instance(j).value(k);
//					distanceSquare += Math.pow((valueA - valueB), 2);
//				}// Of for k
//				distanceMatrix[i][j] = Math.sqrt(distanceSquare);
//			}// Of for j
//		}// Of for i
//	}// Of computeEuclideanDistance
//	
//	public void computeManhattanDistance() {
//		distanceMatrix = new double[numInstances()][numInstances()];
//		for (int i = 0; i < numInstances(); i ++) {
//			for (int j = 0; j < numInstances(); j ++) {
//				double distance = 0;
//				for (int k = 0; k < numAttributes() - 1; k ++) {
//					double valueA = instance(i).value(k);
//					double valueB = instance(j).value(k);
//					distance += Math.abs(valueA - valueB);
//				}// Of for k
//				distanceMatrix[i][j] = distance;
//			}// Of for j
//		}// Of for i
//	}// Of computeManhattanDistance
//	
//	public void computeCosineDistance() {
//		distanceMatrix = new double[numInstances()][numInstances()];
//		for (int i = 0; i < numInstances(); i ++) {
//			for (int j = 0; j < numInstances(); j ++) {
//				double multiResult = 0;
//				double XSquare = 0;
//				double YSquare = 0;
//				for (int k = 0; k < numAttributes() - 1; k ++) {
//					double valueA = instance(i).value(k);
//					double valueB = instance(j).value(k);
//					multiResult += valueA * valueB;
//					XSquare += Math.pow(valueA, 2);
//					YSquare += Math.pow(valueB, 2);
//				}// Of for k
//				distanceMatrix[i][j] = multiResult / (Math.sqrt(XSquare) * Math.sqrt(YSquare));
//			}// Of for j
//		}// Of for i
//	}// Of computeCosineDistance
//	
//	public void getDistanceMatrix() {
//		switch (dataSetType) {
//		case 1:
//			computeBoolDistance();
//			break;
//		case 2:
//			computeEuclideanDistance();
//			break;
//		case 3:
//			computeManhattanDistance();
//			break;
//		case 4:
//			computeCosineDistance();
//			break;
//		default :
//			System.out.println("Error occured in getDistanceMatrix!");
//			break;
//		}// Of switch
//	}// Of getDistanceMatrix
//	
//	public void getMisCostMatrix() {
//		misCostMatrix = new double[numClasses()][numClasses()];
//		misCostMatrix[0][0] = 0;
//		misCostMatrix[1][1] = 0;
//		misCostMatrix[0][1] = misA;
//		misCostMatrix[1][0] = misB;
//	}// Of getMisCostMatrix
//	
//	public void getInsertAttribute() {
//		int classIndex = numAttributes() - 1;
//		Enumeration classValuesEnumeration = attribute(classIndex).enumerateValues();
//		FastVector values = new FastVector();
//		while(classValuesEnumeration.hasMoreElements()) {
//			values.addElement("" + classValuesEnumeration.nextElement());
//		}// Of while
//		String title = "Predict";
//		Attribute insertAttribute = new Attribute(title, values);
//		insertAttributeAt(insertAttribute, 0);
//	}// Of getInsertAttribute
//	
//	public void getAnnounceLength() {
//		instanceTypes = new int[numInstances()];
//		for (int i = 0; i < numInstances(); i ++) {
//			instanceTypes[i] = DELAYED; 
//		}//Of for i
//		
//		expectedCosts = new double[numInstances()];
//	}// Of getAnnounceLength
//	
//	public void getInitialTrainingSet(double paraPercentage) {
//		boolean[] isInitialBought = null;
//		try {
//			isInitialBought = SimpleTool.generateBooleanArrayForDivision(numInstances(), paraPercentage);
//		} catch (Exception e) {
//			System.out.println("Error occured in initializeTaughtInstance!");
//		}// Of try
//		
//		for (int i = 0; i < isInitialBought.length; i ++) {
//			if (isInitialBought[i]) {
//				instance(i).setValue(0, instance(i).value(numAttributes() - 1));
//				instanceTypes[i] = TAUGHT;
//			}// Of if
//		}// Of for i
//	}// Of getInitialTrainingSet
//	
//	public void getInitialization(double paraPercentage) {
//		getDistanceMatrix();
//		getMisCostMatrix();
//		getInsertAttribute();
//		getAnnounceLength();
//		getInitialTrainingSet(paraPercentage);
//	}// Of initialization
//	
//	public void getParameterSetting() {
//		teachCost = 10;
//		misA = 40;
//		misB = 60;
//		dataSetType = 2;
//		kValue = 10;
//		stepLength = 0.02;
//		
//		boughtValue = (int)(stepLength * numInstances());
//		
//		totalTaught = 0;
//		totalGuessed = 0;
//	}// Of getParameterSetting
//	
//	public int[] getNearestNeighbors(int paraIndex) {
//		int[] tempNeighbors = new int[kValue + 1];
//		double[] tempDistances = new double[kValue + 1];
//		int tempDeterminedNeighors = 0;
//		for (int i = 0; i < numInstances(); i ++) {
//			if (instanceTypes[i] == DELAYED) {
//				continue;
//			}//Of if
//			
//			double tempCurrentDistance = distanceMatrix[paraIndex][i];
//			int j;
//			
//			for (j = tempDeterminedNeighors - 1; j >= 0; j --) {
//				if (tempDistances[j] > tempCurrentDistance) {
//					tempDistances[j + 1] = tempDistances[j];
//					tempNeighbors[j + 1] = tempNeighbors[j];
//				} else {
//					break;
//				}//Of if
//			}//Of for j
//			
//			//Insert now
//			tempDistances[j + 1] = tempCurrentDistance;
//			tempNeighbors[j + 1] = i;
//			
//			if (tempDeterminedNeighors < kValue) {
//				tempDeterminedNeighors ++;
//			}//Of if
//		}//Of for i
//		
//		int[] resultNeighbors =  new int[kValue];
//		for (int i = 0; i < kValue; i ++) {
//			resultNeighbors[i] = tempNeighbors[i];
//		}// Of for i
//		
//		return resultNeighbors;
//	}//Of getNearestNeighbors
//	
//	public double getExpectedCost(int[] paraNeighbors) {
//		double tempCost = 0;
//		double tempMinimalCost = 10000;
//		int tempClassIndex = numAttributes() - 1;
//		int tempClassLabels = attribute(tempClassIndex).numValues();
//		int tempNeighborClass;
//		
//		for (int i = 0; i < tempClassLabels; i ++) {
//			tempCost = 0;
//			for (int j = 0; j < paraNeighbors.length; j ++) {
//				tempNeighborClass = (int)instance(paraNeighbors[j]).value(tempClassIndex);
//				tempCost += misCostMatrix[tempNeighborClass][i];
//			}//Of for j
//			
//			if (tempMinimalCost > tempCost) {
//				tempMinimalCost = tempCost;
//				currentBestLabel = i;
//			}//Of if
//		}//Of for i
//		
//		return tempMinimalCost / kValue;
//	}//Of getExpectedCost
//
//	public void getOneRoundGuessedInstances() {
//		int[] tempNeighborIndices;
//		double tempMisCost;
//		currentGuessed = 0;
//
//		for (int i = 0; i < numInstances(); i ++) {
//			if (instanceTypes[i] == TAUGHT) {
//				continue;
//			} else if (instanceTypes[i] == GUESSED){
//				continue;
//			}// Of if
//			
//			tempNeighborIndices = getNearestNeighbors(i);
//			tempMisCost = getExpectedCost(tempNeighborIndices);
//			expectedCosts[i] = tempMisCost;
//			if (tempMisCost < teachCost) {
//				instance(i).setValue(0, currentBestLabel);
//				instanceTypes[i] = GUESSED;
//				currentGuessed ++;
//				totalGuessed ++;
//			}//Of if
//		}//Of for I
//	}//Of getOneRoundGuessedInstances
//	
//	public int getOneRoundTaughtInstances() {
//		int[] tempIndices = new int[boughtValue + 1];
//		double[] tempCosts = new double[boughtValue + 1];
//		currentTaught = 0;
//		
//		for (int i = 0; i < numInstances(); i ++) {
//			if (instanceTypes[i] != DELAYED) {
//				continue;
//			}//Of if
//		
//			double tempExceptedMisclassificationCosts = expectedCosts[i];
//			int j;
//			
//			for (j = currentTaught - 1; j >= 0; j--) {
//				if (tempCosts[j] < tempExceptedMisclassificationCosts) {
//					tempCosts[j + 1] = tempCosts[j];
//					tempIndices[j + 1] = tempIndices[j];
//				} else {
//					break;
//				}// Of if 
//			}// Of for j
//			
//			tempCosts[j + 1] = tempExceptedMisclassificationCosts;
//			tempIndices[j + 1] = i;
//			
//			if (currentTaught < boughtValue) {
//				currentTaught ++;
//			}// Of if 
//		}//Of for i
//		
//		for (int i = 0; i < currentTaught; i ++) {
//			instanceTypes[tempIndices[i]] = TAUGHT;
//			//Copy the actual class label.
//			instance(tempIndices[i]).setValue(0, instance(tempIndices[i]).value(numAttributes() - 1));
//			totalTaught ++;
//		}//Of for i
//		
//		return currentTaught;
//	}//Of getOneRoundTaughtInstances
//	
//	public int[][] getIterativeLearning() {
//		int[][] tempRecordMatrix = new int[numInstances() / boughtValue][2];
//		int tempTimes = 0;
//		
//		boolean isTerminal = getTerminationSign();
//		while (!isTerminal) {
//			getOneRoundGuessedInstances();
//			tempRecordMatrix[tempTimes][0] = currentGuessed;
//			getOneRoundTaughtInstances();
//			tempRecordMatrix[tempTimes][1] = currentTaught;
//			isTerminal = getTerminationSign();
//			tempTimes ++;
//		}// Of while
//		
//		int[][] recordMatrix = new int[tempTimes][2];
//		for (int i = 0; i < tempTimes; i ++) {
//			recordMatrix[i] = tempRecordMatrix[i];
//		}// Of for i
//		
//		return recordMatrix;
//	}// Of getIterativeLearning
//	
//	public boolean getTerminationSign() {
//		int initialInstance = (int)(percentage * numInstances());
//		int total = initialInstance + totalGuessed + totalTaught;
//		
//		boolean flag = true;
//		if(total != numInstances()) {
//			flag = false;
//		}// Of 
//		
//		return flag;
//	}// Of getTerminationSign
//	
//	public int[] getVerifiedResult() {
//		int tempCorrect = 0;
//		int tempA = 0;
//		int tempB = 0;
//		int[] tempClassificationResult = new int[3];
//		for (int i = 0; i < numInstances(); i ++) {
//			if (instanceTypes[i] == GUESSED) {
//				int tempPredictValue = (int)instance(i).value(0);
//				int tempTrueValue = (int)instance(i).value(numAttributes() - 1);
//				if (tempPredictValue == tempTrueValue) {
//					tempCorrect ++;
//				} else {
//					if (tempTrueValue == 0) {
//						tempA ++;
//					} else {
//						tempB ++;
//					}// Of if
//				}// Of if
//			}// Of if
//		}// Of for i
//		
//		tempClassificationResult[0] = tempCorrect;
//		tempClassificationResult[1] = tempA;
//		tempClassificationResult[2] = tempB;
//		
//		return tempClassificationResult;
//	}// Of getVerifiedResult
//	
//	public static int[][][] learnerTest() {
//		Learner tempLearner = null;
//		Learner firstSet = null;
//		Learner secondSet = null;
//		Learner thirdSet = null;
//		
//		String arffFilename = "E:/data/numerical/nominalize/EE.arff";
//		
//		try {
//			FileReader fileReader = new FileReader(arffFilename);
//			tempLearner = new Learner(fileReader);
//			fileReader.close();
//			
//			firstSet = new Learner(tempLearner);
//			secondSet = new Learner(tempLearner);
//			thirdSet = new Learner(tempLearner);
//		} catch (Exception ee) {
//			System.out.println("Error occured in learnerTest!");
//		}// Of try
//		
//		int[][] noneResult = compareNone(firstSet);
//		int[][] kNNResult = compareKNN(secondSet);	
//		int[][] TreeResult = compareDecisionTree(thirdSet);
//		
//		int[][][] printResult = new int[3][][];
//		printResult[0] = kNNResult;
//		printResult[1] = TreeResult;
//		printResult[2] = noneResult;
//		
//		return printResult;
//	}// Of labelTest
//
//	public static void main(String[] args) {
//		try {
//			FileOutputStream fs = new FileOutputStream(new File("D:\\EE10.txt"));
//			PrintStream p = new PrintStream(fs);
//			
//			for (int i = 0; i < 10; i ++) {
//				percentage = (i + 1) * 0.01;
//				p.println(percentage);
//				for (int j = 0; j < 10; j ++) {
//					int[][][] tempCubic = learnerTest();
//					p.println(Arrays.deepToString(tempCubic));
//				}// Of for j
//			}// Of for i
//			
//			p.close();
//		} catch (Exception e) {
//			System.out.println("Error occured in main!");
//			e.printStackTrace();
//		}// Of try
//		
//		System.out.println("Finish!");
//	}// Of main
//
//	public static int[][] compareNone(Learner tempLearner) {
//		tempLearner.getParameterSetting();
//		tempLearner.getInitialization(percentage);
//		int[][] recordMatrix = tempLearner.getIterativeLearning();
//		int[] correctAndUncorrect = tempLearner.getVerifiedResult();
//		
//		int[][] returnResult = new int[recordMatrix.length + 2][2];
//		returnResult[0] = correctAndUncorrect;
//		returnResult[1][0] = tempLearner.totalGuessed;
//		returnResult[1][1] = tempLearner.totalTaught;
//		for (int i = 2; i < returnResult.length; i ++) {
//			returnResult[i] = recordMatrix[i - 2];
//		}// Of for i
//		
//		noneTotalTaught = tempLearner.totalTaught;
//		
//		return returnResult;
//	}// Of compareSelf
//	
//	public int getNoneClassLabel(int paraIndex) {
//		double[] distanceVector = distanceMatrix[paraIndex];
//		
//		int count = 0;
//		for (int i = 0; i < numInstances(); i ++) {
//			if (instanceTypes[i] == TAUGHT) {
//				count ++;
//			}// Of if
//		}// Of for i
//		
//		int index = 0;
//		double[][] usefulInformation = new double[count][2];
//		for (int i = 0; i < numInstances(); i ++) {
//			if (instanceTypes[i] == TAUGHT) {
//				usefulInformation[index][0] = i;
//				usefulInformation[index][1] = distanceVector[i];
//				index ++;
//			}// Of if
//		}// Of for i
//		
//		for (int i = 0; i < usefulInformation.length; i ++) {
//			for (int j = i + 1; j < usefulInformation.length; j ++) {
//				if (usefulInformation[i][1] > usefulInformation[j][1]) {
//					double[] tempArray = usefulInformation[i];
//					usefulInformation[i] = usefulInformation[j];
//					usefulInformation[j] = tempArray;
//				}// Of if
//			}// Of for j
//		}// Of for i
//		
//		int[] neighborIndices = new int[kValue];
//		int[] vote = new int[numClasses()];
//		for (int i = 0; i < neighborIndices.length; i ++) {
//			int classLabel = (int)instance(neighborIndices[i]).classValue();
//			vote[classLabel] ++;
//		}// Of for i
//		
//		int max = Integer.MIN_VALUE;
//		int predictedLabel = -1;
//		for (int i = 0; i < vote.length; i ++) {
//			if (vote[i] > max) {
//				max = vote[i];
//				predictedLabel = i;
//			}// Of if
//		}// Of for i
//		
//		return predictedLabel;
//	}// Of getNoneClassLabel
//	
//	public void getKNNClassification() {
//		for (int i = 0; i < numInstances(); i ++) {
//			if (instanceTypes[i] == TAUGHT) {
//				continue;
//			}// Of if
//			
//			int predictedLabel = getNoneClassLabel(i);
//			instance(i).setValue(0, predictedLabel);
//		}// Of for i
//	}// Of getKNNClassification
//	
//	public int[] getNoneVerifiedResult() {
//		int tempCorrect = 0;
//		int tempA = 0;
//		int tempB = 0;
//		
//		for (int i = 0; i < numInstances(); i ++) {
//			int predictedLabel = (int)instance(i).value(0);
//			int intrinsicLabel = (int)instance(i).classValue();
//			if (predictedLabel == intrinsicLabel) {
//				tempCorrect ++;
//			} else if (predictedLabel != intrinsicLabel) {
//				if (intrinsicLabel == 0) {
//					tempA ++;
//				} else {
//					tempB ++;
//				}// Of if
//			} else {
//				System.out.println("Error occured in getNoneVerifiedResult!");
//			}// Of if
//		}// Of for i
//
//		int[] tempClassificationResult = new int[3];
//		tempClassificationResult[0] = tempCorrect;
//		tempClassificationResult[1] = tempA;
//		tempClassificationResult[2] = tempB;
//		
//		return tempClassificationResult;
//	}// Of getNoneVerifiedResult
//	
//	public static int[][] compareKNN(Learner tempLearner) {
//		int labeled = (int)(percentage * tempLearner.numInstances()) + noneTotalTaught;
//		double newPercentage = (labeled + 0.0) / tempLearner.numInstances();
//		
//		tempLearner.getParameterSetting();
//		tempLearner.getInitialization(newPercentage);
//		tempLearner.getKNNClassification();
//		int[] classificationResult = tempLearner.getNoneVerifiedResult();
//		int[][] returnResult = new int[1][];
//		returnResult[0] = classificationResult;
//		
//		return returnResult;
//	}// Of compareKNN
//	
//	public Learner[] getDividedDataSets(double paraPercentage) throws Exception {
//		boolean[] firstInclusionArray = SimpleTool.generateBooleanArrayForDivision(numInstances(), paraPercentage);
//		Learner firstDecisionSystem = new Learner(this);
//		firstDecisionSystem.delete(firstInclusionArray);
//
//		boolean[] secondInclusionArray = SimpleTool.revertBooleanArray(firstInclusionArray);
//		Learner secondDecisionSystem = new Learner(this);
//		secondDecisionSystem.delete(secondInclusionArray);
//
//		Learner[] subsets = new Learner[2];
//		subsets[0] = firstDecisionSystem;
//		subsets[1] = secondDecisionSystem;
//		
//		return subsets;
//	}// Of getDividedDataSets
//	
//	public static int[][] compareDecisionTree(Learner tempLearner) {
//		int labeled = (int)(percentage * tempLearner.numInstances()) + noneTotalTaught;
//		double newPercentage = (labeled + 0.0) / tempLearner.numInstances();
//		
//		Learner trainingSet = null;
//		Learner testingSet = null;
//		Learner[] tempSets = null;
//		
//		Learner trainingCSId3 = null;
//		Learner trainingCSJ48 = null;
//		Learner testingCSId3 = null;
//		Learner testingCSJ48 = null;
//		
//		try {
//			tempSets = tempLearner.getDividedDataSets(newPercentage);
//			trainingSet = tempSets[0];
//			testingSet = tempSets[1];
//			
//			trainingCSId3 = new Learner(trainingSet);
//			trainingCSJ48 = new Learner(trainingSet);
//			testingCSId3 = new Learner(testingSet);
//			testingCSJ48 = new Learner(testingSet);
//		} catch (Exception e) {
//			e.printStackTrace();
//			System.out.println("Error occured in compareDecisionTree!");
//		}// Of try
//		
//		//int[] CSId3Result = compareCSId3(trainingCSId3, testingCSId3);
//		int[] CSJ48Result = compareCSJ48(trainingCSJ48, testingCSJ48);
//		//int[][] returnResult = new int[2][];
//		//returnResult[0] = CSId3Result;
//		//returnResult[1] = CSJ48Result;
//		int[][] returnResult = new int[2][];
//		//returnResult[0] = CSId3Result;
//		returnResult[0] = null;
//		returnResult[1] = CSJ48Result;
//		
//		return returnResult;
//	}// Of compareDecisionTrees
//	
//	public Learner getCSId3TrainingSet() {
//		Learner CSTrainingSet = null;
//		try {
//			CSTrainingSet = new Learner(this);
//		} catch (Exception e) {
//			e.printStackTrace();
//			System.out.println("Error occured in getCSId3TrainingSet!");
//		}// Of try
//		
//		CSTrainingSet.deleteWithMissingClass();
//		Id3 classifier = new Id3();
//		
//		costMatrixWeka = new CostMatrix(CSTrainingSet.numClasses());
//		costMatrixWeka.setCell(0, 0, 0);
//		costMatrixWeka.setCell(0, 1, misA);
//		costMatrixWeka.setCell(1, 0, misB);
//		costMatrixWeka.setCell(1, 1, 0);
//	
//		Random random = null;
//		if (!(classifier instanceof WeightedInstancesHandler)) {
//			random = new Random(1);
//		}// Of if
//		try {
//			Instances CSDataSet = costMatrixWeka.applyCostMatrix(CSTrainingSet, random);
//			CSTrainingSet = new Learner(CSDataSet);
//		} catch (Exception e) {
//			e.printStackTrace();
//			System.out.println("Error occured in getCSId3TrainingSet!");
//		}// Of try
//		
//		return CSTrainingSet;
//	}// Of getCSId3TrainingSet
//	
//	public static int[] compareCSId3(Learner trainingSet, Learner testingSet) {
//		Learner CSTrainingSet = trainingSet.getCSId3TrainingSet();
//		
//		int tempCorrect = 0;
//		int tempA = 0;
//		int tempB = 0;
//		int tempUnclassified = 0;
//		
//		Id3 classifier = new Id3();
//		try {
//			classifier.buildClassifier(CSTrainingSet);
//			for (int i = 0; i < testingSet.numInstances(); i ++) {
//				int predictedLabel = (int)classifier.classifyInstance(testingSet.instance(i));
//				int trueLabel = (int)testingSet.instance(i).classValue();
//				if (trueLabel < 0 || trueLabel >= testingSet.numClasses()) {
//					tempUnclassified ++;
//				} else {
//					if (trueLabel == predictedLabel) {
//						tempCorrect ++;
//					} else {
//						if (trueLabel == 0) {
//							tempA ++;
//						} else {
//							tempB ++;
//						}// Of if
//					}// Of if
//				}// Of if
//			}// Of for i
//		} catch (Exception e) {
//			e.printStackTrace();
//			System.out.println("Error occured in compareCSId3!");
//		}// Of try
//		
//		int[] returnResult = new int[4];
//		returnResult[0] = tempCorrect;
//		returnResult[1] = tempA;
//		returnResult[2] = tempB;
//		returnResult[3] = tempUnclassified;
//		
//		return returnResult;
//	}// Of compareCSId3
//	
//	public Learner getCSJ48TrainingSet() {
//		Learner CSTrainingSet = null;
//		try {
//			CSTrainingSet = new Learner(this);
//		} catch (Exception e) {
//			e.printStackTrace();
//			System.out.println("Error occured in getCSJ48TrainingSet!");
//		}// Of try
//		
//		CSTrainingSet.deleteWithMissingClass();
//		J48 classifier = new J48();
//		
//		costMatrixWeka = new CostMatrix(CSTrainingSet.numClasses());
//		costMatrixWeka.setCell(0, 0, 0);
//		costMatrixWeka.setCell(0, 1, misA);
//		costMatrixWeka.setCell(1, 0, misB);
//		costMatrixWeka.setCell(1, 1, 0);
//	
//		Random random = null;
//		if (!(classifier instanceof WeightedInstancesHandler)) {
//			random = new Random(1);
//		}// Of if
//		try {
//			Instances CSDataSet = costMatrixWeka.applyCostMatrix(CSTrainingSet, random);
//			CSTrainingSet = new Learner(CSDataSet);
//		} catch (Exception e) {
//			e.printStackTrace();
//			System.out.println("Error occured in getCSId3TrainingSet!");
//		}// Of try
//		
//		return CSTrainingSet;
//	}// Of getCSJ48TrainingSet
//	
//	public static int[] compareCSJ48(Learner trainingSet, Learner testingSet) {
//		int tempCorrect = 0;
//		int tempA = 0;
//		int tempB = 0;
//		int tempUnclassified = 0;
//		
//		J48 classifier = new J48();
//		try {
//			classifier.buildClassifier(trainingSet);
//			for (int i = 0; i < testingSet.numInstances(); i ++) {
//				int predictedLabel = (int)classifier.classifyInstance(testingSet.instance(i));
//				int trueLabel = (int)testingSet.instance(i).classValue();
//				if (trueLabel < 0 || trueLabel >= testingSet.numClasses()) {
//					tempUnclassified ++;
//				} else {
//					if (trueLabel == predictedLabel) {
//						tempCorrect ++;
//					} else {
//						if (trueLabel == 0) {
//							tempA ++;
//						} else {
//							tempB ++;
//						}// Of if
//					}// Of if
//				}// Of if
//			}// Of for i
//		} catch (Exception e) {
//			e.printStackTrace();
//			System.out.println("Error occured in compareCSJ48!");
//		}// Of try
//		
//		int[] returnResult = new int[4];
//		returnResult[0] = tempCorrect;
//		returnResult[1] = tempA;
//		returnResult[2] = tempB;
//		returnResult[3] = tempUnclassified;
//		
//		return returnResult;
//	}// Of compareCSJ48
//}// Of class Learner



import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Reader;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Random;

import coser.common.SimpleTool;
import coser.datamodel.decisionsystem.CostSensitiveDecisionSystem;
import weka.classifiers.CostMatrix;
import weka.classifiers.trees.Id3;
import weka.classifiers.trees.J48;
import weka.core.*;

public class Learner extends CostSensitiveDecisionSystem {
	
	private static final long serialVersionUID = -82097437233430226L;
	
	public static double percentage;
	public static double stepLength;
	public static int noneTotalTaught;
	
	public static double teachCost;
	public static double misCostA;
	public static double misCostB;
	public int dataSetType;
	public int kValue;
	public int boughtValue;
	
	private int currentBestLabel;
	private int totalGuessed;
	private int totalTaught;
	private int currentGuessed;
	private int currentTaught;
	
	public double[][] misCostMatrix;
	public int[] instanceTypes;
	public double[] expectedCosts;
	public double[][] distanceMatrix;
	
	public CostMatrix costMatrixWeka;
		
	public static final int TAUGHT  = 0;
	public static final int GUESSED = 1;
	public static final int DELAYED = 2;
	
	
	
	/*****
	 * Constructor One
	 * @param paraInstances
	 * @throws IOException
	 * @throws Exception
	 */
	public Learner(Instances paraInstances) throws IOException, Exception {
		super(paraInstances);
	}// Of the first constructor

	
	
	/*****
	 * Constructor Two
	 * @param paraReader
	 * @throws IOException
	 * @throws Exception
	 */
	public Learner(Reader paraReader) throws IOException, Exception {
		super(paraReader);
	}// Of the second constructor

	
	
	/*****
	 * Compute the Overlap distance between each pair of instances.distance越小越相似
	 * For symbolic data.
	 */
	public void computeOverlapDistance() {
		// initial member variable [distanceMatrix]
		distanceMatrix = new double[numInstances()][numInstances()];
		
		for (int i = 0; i < numInstances(); i ++) {
			for (int j = 0; j < numInstances(); j ++) {
				// temporary variable
				double distance = 0;
				for (int k = 0; k < numAttributes() - 1; k ++) {
					// get the value of attribute
					double valueA = instance(i).value(k);
					double valueB = instance(j).value(k);
					// match the values
					if ((int)valueA != (int)valueB) {
						distance ++;
					}// Of if
				}// Of for k
				// assignment
				distanceMatrix[i][j] = distance;
			}// Of for j
		}// Of for i
	}// Of computeOverlapDistance
	
	
	
	/*****
	 * Compute the Euclidean distance between each pair of instances.
	 * For numeric data.
	 */
	public void computeEuclideanDistance() {
		// initial member variable [distanceMatrix]
		distanceMatrix = new double[numInstances()][numInstances()];
		
		for (int i = 0; i < numInstances(); i ++) {
			for (int j = 0; j < numInstances(); j ++) {
				// temporary variable
				double distanceSquare = 0;
				for (int k = 0; k < numAttributes() - 1; k ++) {
					// get the value of attribute
					double valueA = instance(i).value(k);
					double valueB = instance(j).value(k);
					distanceSquare += Math.pow((valueA - valueB), 2);
				}// Of for k
				// assignment
				distanceMatrix[i][j] = Math.sqrt(distanceSquare);
			}// Of for j
		}// Of for i
	}// Of computeEuclideanDistance
	
	
	
	/*****
	 * Compute the Manhattan distance between each pair of instances.
	 * For numeric data.
	 */
	public void computeManhattanDistance() {
		// initial member variable [distanceMatrix]
		distanceMatrix = new double[numInstances()][numInstances()];
		
		for (int i = 0; i < numInstances(); i ++) {
			for (int j = 0; j < numInstances(); j ++) {
				// temporary variable
				double distance = 0;
				for (int k = 0; k < numAttributes() - 1; k ++) {
					// get the value of attribute
					double valueA = instance(i).value(k);
					double valueB = instance(j).value(k);
					distance += Math.abs(valueA - valueB);
				}// Of for k
				// assignment
				distanceMatrix[i][j] = distance;
			}// Of for j
		}// Of for i
	}// Of computeManhattanDistance
	
	
	
	/*****
	 * Compute the Cosine distance between each pair of instances.
	 * For numeric data.
	 */
	public void computeCosineDistance() {
		// initial member variable [distanceMatrix]
		distanceMatrix = new double[numInstances()][numInstances()];
		
		for (int i = 0; i < numInstances(); i ++) {
			for (int j = 0; j < numInstances(); j ++) {
				// temporary variables
				double multiResult = 0;
				double XSquare = 0;
				double YSquare = 0;
				for (int k = 0; k < numAttributes() - 1; k ++) {
					// get the value of attribute
					double valueA = instance(i).value(k);
					double valueB = instance(j).value(k);
					multiResult += valueA * valueB;
					XSquare += Math.pow(valueA, 2);
					YSquare += Math.pow(valueB, 2);
				}// Of for k
				// assignment
				distanceMatrix[i][j] = multiResult / (Math.sqrt(XSquare) * Math.sqrt(YSquare));
			}// Of for j
		}// Of for i
	}// Of computeCosineDistance
	
	
	
	/*****
	 * Select the distance measure.
	 */
	public void getDistanceMatrix() {
		switch (dataSetType) {
		case 1:
			computeOverlapDistance();
			break;
		case 2:
			computeEuclideanDistance();
			break;
		case 3:
			computeManhattanDistance();
			break;
		case 4:
			computeCosineDistance();
			break;
		default :
			System.out.println("Error occured in getDistanceMatrix!");
			break;
		}// Of switch
	}// Of getDistanceMatrix
	
	
	
	/*****
	 * Set misclassification cost matrix
	 */
	public void getMisCostMatrix() {
		// initial member variable [misCostMatrix]
		misCostMatrix = new double[numClasses()][numClasses()];
		
		misCostMatrix[0][0] = 0;
		misCostMatrix[0][1] = misCostA;
		misCostMatrix[1][0] = misCostB;
		misCostMatrix[1][1] = 0;
		
//		for (int i = 0; i < numClasses(); i ++) {
//			for (int j = 0; j < numClasses(); j ++) {
//				if (i == j) {
//					misCostMatrix[i][j] = 0;
//				} else {
//					misCostMatrix[i][j] = misCost;
//				}// Of if
//			}// Of for j
//		}// Of for i
	}// Of getMisCostMatrix
	
	
	
	/*****
	 * Insert a new attribute to mark statue
	 */
	public void getInsertAttribute() {
		// locate the decision attribute
		int classIndex = numAttributes() - 1;
		// copy the domain of decision attribute
		Enumeration classValuesEnumeration = attribute(classIndex).enumerateValues();
		FastVector values = new FastVector();
		// assignment
		while(classValuesEnumeration.hasMoreElements()) {
			values.addElement("" + classValuesEnumeration.nextElement());
		}// Of while
		// rename
		String title = "Predict";
		// create object
		Attribute insertAttribute = new Attribute(title, values);
		insertAttributeAt(insertAttribute, 0);
	}// Of getInsertAttribute
	
	
	
	/*****
	 * Initial the statue of instances.
	 */
	public void initalStatue() {
		// initial member variable [instanceType]
		instanceTypes = new int[numInstances()];
		
		for (int i = 0; i < numInstances(); i ++) {
			instanceTypes[i] = DELAYED; 
		}//Of for i
		
		// initial member variable [expectedCosts]
		expectedCosts = new double[numInstances()];
	}// Of initalStatue
	
	
	
	/*****
	 * Generate training set.
	 * @param paraPercentage the scale of training set.
	 */
	public void getInitialTrainingSet(double paraPercentage) {
		// temporary boolean vector for recording the labeled instances
		boolean[] isInitialBought = null;
		
		try {
			// generate random vector
			isInitialBought = SimpleTool.generateBooleanArrayForDivision(numInstances(), paraPercentage);
		} catch (Exception e) {
			System.out.println("Error occured in initializeTaughtInstance!");
		}// Of try
		
		for (int i = 0; i < isInitialBought.length; i ++) {
			if (isInitialBought[i]) {
				// assignment on new added attribute？？？？
				instance(i).setValue(0, instance(i).value(numAttributes() - 1));
				// change the statue of target instance
				instanceTypes[i] = TAUGHT;
			}// Of if
		}// Of for i
	}// Of getInitialTrainingSet
	
	
	
	/*****
	 * Get the initial information.
	 * @param paraPercentage the scale of training set.
	 */
	public void getInitialization(double paraPercentage) {
		getDistanceMatrix();
		getMisCostMatrix();
		getInsertAttribute();
		initalStatue();
		getInitialTrainingSet(paraPercentage);
	}// Of initialization
	
	
	
	/*****
	 * Set the initial parameters.
	 */
	public void getParameterSetting() {
		teachCost = 10;     // the setting of teacher cost
		
//==========================================		
//		misCostA = 90;
//		misCostB = 10;
//==========================================
		
		dataSetType = 1;    // the type of distance measure
		kValue = 10;        // the setting of k
		stepLength = 0.02;  // the setting of step
		boughtValue = (int)(stepLength * numInstances());//boughtValue=8;numInstances()=435
		System.out.println("boughtValue"+boughtValue);
		
		totalTaught = 0;    // the sum number of taught instances
		totalGuessed = 0;   // the sum number of guessed instances
	}// Of getParameterSetting
	
	
	
	/*****
	 * Get the nearest neighbors of this unlabeled instance.
	 * @param paraIndex the current unlabeled instance.
	 * @return resultNeighbors the vector of nearest neighbors.
	 */
	public int[] getNearestNeighbors(int paraIndex) {
		// temporary integer vector for recording the indices of neighbors
		int[] tempNeighbors = new int[kValue + 1];
		// temporary double vector for recording the distances of neighbors
		double[] tempDistances = new double[kValue + 1];
		
		int tempDeterminedNeighors = 0;
		for (int i = 0; i < numInstances(); i ++) {
			// ignore other unlabeled instances，计算未分类实例与已分类实例间的距离
			if (instanceTypes[i] == DELAYED) {
				continue;
			}//Of if
			//距离矩阵在初始化时已经计算完成
			double tempCurrentDistance = distanceMatrix[paraIndex][i];
			int j;			
			//从后往前比较
			for (j = tempDeterminedNeighors - 1; j >= 0; j --) {
				if (tempDistances[j] > tempCurrentDistance) {
					tempDistances[j + 1] = tempDistances[j];
					tempNeighbors[j + 1] = tempNeighbors[j];
				} else {
					break;
				}//Of if
			}//Of for j
			//Insert now,按距离从小到大排列
			tempDistances[j + 1] = tempCurrentDistance;
			tempNeighbors[j + 1] = i;
			
			if (tempDeterminedNeighors < kValue) {
				tempDeterminedNeighors ++;
			}//Of if
		}//Of for i
		
		int[] resultNeighbors =  new int[kValue];
		for (int i = 0; i < kValue; i ++) {
			resultNeighbors[i] = tempNeighbors[i];
		}// Of for i
		
		return resultNeighbors;
	}//Of getNearestNeighbors
	
	public double getExpectedCost(int[] paraNeighbors) {
		double tempCost = 0;
		double tempMinimalCost = 10000;
		int tempClassIndex = numAttributes() - 1;
		int tempClassLabels = attribute(tempClassIndex).numValues();
		
		int tempNeighborClass;
		
		for (int i = 0; i < tempClassLabels; i ++) {
			tempCost = 0;
			//根据邻居的label与当前预测的label，结合代价矩阵，算出一个cost
			for (int j = 0; j < paraNeighbors.length; j ++) {
				tempNeighborClass = (int)instance(paraNeighbors[j]).value(tempClassIndex);
				tempCost += misCostMatrix[tempNeighborClass][i];
			}//Of for j
			//期望值越小越好
			if (tempMinimalCost > tempCost) {
				tempMinimalCost = tempCost;
				currentBestLabel = i;
			}//Of if
		}//Of for i
		
		return tempMinimalCost / kValue;
	}//Of getExpectedCost

	
	
	/*****
	 * Find the instances which can be guessed in this iteration.
	 */
	public void getOneRoundGuessedInstances() {
		// temporary integer vector for recording the neighbors
		int[] tempNeighborIndices;
		// temporary double for recording the misclassification cost of current instance
		double tempMisCost;
		// the number of guessed instances in this iteration
		currentGuessed = 0;

		for (int i = 0; i < numInstances(); i ++) {
			// Step 1. ignore labeled instances
			if (instanceTypes[i] == TAUGHT) {
				continue;
			} else if (instanceTypes[i] == GUESSED){
				continue;
			}// Of if
			
			// Step 2. get nearest neighbors of this unlabeled instance
			tempNeighborIndices = getNearestNeighbors(i);
			
			// Step 3. compute the expected misclassification cost of this unlabeled instance
			tempMisCost = getExpectedCost(tempNeighborIndices);
			// Step 3.1 record the expected misclassification cost
			expectedCosts[i] = tempMisCost;
			
			// Step 4. measure this unlabeled instance
			if (tempMisCost < teachCost) {
				// set the value of new added attribute
				instance(i).setValue(0, currentBestLabel);
				// change the statue of this unlabeled instance
				instanceTypes[i] = GUESSED;
				
				currentGuessed ++;
				totalGuessed ++;
			}//Of if
		}//Of for I
	}//Of getOneRoundGuessedInstances
	
	
	
	/*****
	 * Find the instances which can be taught in this iteration.
	 * @return currentTaught the number of taught instances in this iteration.
	 */
	public int getOneRoundTaughtInstances() {
		// temporary integer vector for recording indices
		int[] tempIndices = new int[boughtValue + 1];
		// temporary double vector for recording expected misclassification cost
		double[] tempCosts = new double[boughtValue + 1];
		currentTaught = 0;
		
		for (int i = 0; i < numInstances(); i ++) {
			if (instanceTypes[i] != DELAYED) {
				continue;
			}//Of if
			
			//在每一轮迭代中，都会根据kValue个邻居的分类情况重新计算代价期望，
		
			double tempExceptedMisclassificationCosts = expectedCosts[i];
			int j;
			//按照期望从大到小排列
			for (j = currentTaught - 1; j >= 0; j--) {
				if (tempCosts[j] < tempExceptedMisclassificationCosts) {
					tempCosts[j + 1] = tempCosts[j];
					tempIndices[j + 1] = tempIndices[j];
				} else {
					break;
				}// Of if 
			}// Of for j
			//insert
			tempCosts[j + 1] = tempExceptedMisclassificationCosts;
			tempIndices[j + 1] = i;
			
			if (currentTaught < boughtValue) {
				currentTaught ++;
			}// Of if 
		}//Of for i
		
		for (int i = 0; i < currentTaught; i ++) {
			// change the statue of this unlabeled instance
			instanceTypes[tempIndices[i]] = TAUGHT;
			//Copy the actual class label.
			instance(tempIndices[i]).setValue(0, instance(tempIndices[i]).value(numAttributes() - 1));
			totalTaught ++;
		}//Of for i
		
		return currentTaught;
	}//Of getOneRoundTaughtInstances
	
	
	
	/*****
	 * Iterative active learning.
	 * @return recordMatrix
	 */
	public int[][] getIterativeLearning() {
		// temporary integer matrix to record the result on each iteration
		//boughtValue
		int[][] tempRecordMatrix = new int[numInstances() / boughtValue][2];//tempRecordMatrix.length = 54	
		// temporary integer to record the time of iteration
		int tempTimes = 0;
		
		// get the termination of iteration
		boolean isTerminal = getTerminationSign();
		
		// iteration
		while (!isTerminal) {
			// select the instances which can be guessed
			getOneRoundGuessedInstances();
			// record the number of guessed instance in this iteration 
			tempRecordMatrix[tempTimes][0] = currentGuessed;
			
			// select the instances which can be taught
			getOneRoundTaughtInstances();
			// record the number of guessed instance in this iteration
			tempRecordMatrix[tempTimes][1] = currentTaught;
			
			isTerminal = getTerminationSign();
			tempTimes ++;
		}// Of while
		
		// compress
		int[][] recordMatrix = new int[tempTimes][2];
		for (int i = 0; i < tempTimes; i ++) {
			recordMatrix[i] = tempRecordMatrix[i];
		}// Of for i
		
		return recordMatrix;
	}// Of getIterativeLearning
	
	
	
	/*****
	 * Get the termination signal.
	 * @return flag represent the statue of termination.
	 */
	public boolean getTerminationSign() {
		// get the number of training instances，初始状态是学习查询到initialInstance个标签
		int initialInstance = (int)(percentage * numInstances());
		// count the number of labeled instances
		int total = initialInstance + totalGuessed + totalTaught;
		
		boolean flag = true;
		if(total != numInstances()) {
			flag = false;
		}// Of 
		
		return flag;
	}// Of getTerminationSign
	
	
	
	/*****
	 * Compute the classification results.
	 * @return tempClassificationResult the correct and incorrect statistic data.
	 */
	public int[] getVerifiedResult() {
		int tempCorrect = 0;
		int tempUncorrect = 0;
		int tempUncorrectA = 0;
		int tempUncorrectB = 0;
		
		int[] tempClassificationResult = new int[4];
		for (int i = 0; i < numInstances(); i ++) {
			if (instanceTypes[i] == GUESSED) {
				// get the value of predicted label
				int tempPredictValue = (int)instance(i).value(0);
				// get the value of actual label
				int tempTrueValue = (int)instance(i).value(numAttributes() - 1);
				
				// match
				if (tempPredictValue == tempTrueValue) {
					tempCorrect ++;
				} else {
					if (tempTrueValue == 0) {
						
						tempUncorrectA ++;
					} else if (tempTrueValue == 1) {
						tempUncorrectB ++;
					} else {
						
						System.out.println("Error occured in getVerifiedResult!");
					}// Of if
					
					tempUncorrect ++;
				}// Of if
			}// Of if
		}// Of for i
		
		tempClassificationResult[0] = tempCorrect;
		tempClassificationResult[1] = tempUncorrect;
		tempClassificationResult[2] = tempUncorrectA;
		tempClassificationResult[3] = tempUncorrectB;
		
		return tempClassificationResult;
	}// Of getVerifiedResult
	
	
	
	/*****
	 * Experiment.
	 * @return comparison experiment results.
	 */
	public static int[][][] learnerTest() {
		Learner tempLearner = null;
		Learner firstSet = null;
		Learner secondSet = null;
		Learner thirdSet = null;
		
		//String fileName = "voting";
		String arffAdress = "/Users/dengsiyu/eclipse/workspace/Coser_Triple/data/voting.arff";
		
		try {
			// load data
			FileReader fileReader = new FileReader(arffAdress);
			tempLearner = new Learner(fileReader);
			fileReader.close();
			
			// initial objects for experiments
			firstSet = new Learner(tempLearner);
			secondSet = new Learner(tempLearner);
			thirdSet = new Learner(tempLearner);
		} catch (Exception ee) {
			System.out.println("Error occured in learnerTest!");
		}// Of try
		
		int[][] noneResult = null;
		int[][] kNNResult = null;
		int[][] TreeResult = null;
		
		try {
			noneResult = compareNone(firstSet);
			kNNResult  = compareKNN(secondSet);	
			TreeResult = compareDecisionTree(thirdSet);
		} catch (Exception e) {
			System.out.println("Error occured in LearnerTest!");
			e.printStackTrace();
		}// Of try
		
		// collect experiment data
		int[][][] printResult = new int[3][][];
		printResult[0] = kNNResult;
		printResult[1] = TreeResult;
		printResult[2] = noneResult;
		
		return printResult;
	}// Of labelTest

	
	
	/*****
	 * main function.
	 * @param args.
	 */
	public static void main(String[] args) {
		try {
			//haha
			FileOutputStream fs = new FileOutputStream(new File("/Users/dengsiyu/eclipse/workspace/smaleCode/arff/MushroomDD.txt"));
			PrintStream printStrem = new PrintStream(fs);
			printStrem.println(fs);
			// initial member variables [misCostA] and [misCostB]
			misCostA = 100;
			misCostB = 0;	
			System.out.println("第一个二维数组KNN，第二个二维数组：tree；第三个");
			// change the settings of misclassification
			for (int j = 1; j < 10; j ++) {
				misCostA = 100 - j * 10;//90
				misCostB = j * 10;//10
				
				// change the training set
				for (int i = 0; i < 10; i ++) {
					percentage = (i + 1) * 0.01;
					printStrem.println(percentage);
					
					// experiment
					for (int k = 0; k < 10; k ++) {
						int[][][] tempCubic = learnerTest();
						printStrem.println(Arrays.deepToString(tempCubic));
						System.out.println("one Finish");
					}// Of for k
				}// Of for i
				
			}// Of for j
		
			printStrem.close();
		} catch (Exception e) {
			System.out.println("Error occured in main!");
			e.printStackTrace();
		}// Of try
		
		System.out.println("Finish!");
	}// Of main

	
	
	/*****
	 * Compare with itself.
	 * @param tempLearner the basic decision system.
	 * @return classification results.
	 */
	public static int[][] compareNone(Learner tempLearner) {
		// Step 1. initial settings
		//teachCost Kvalue boughtValue dataSetType
		tempLearner.getParameterSetting();
		//距离矩阵，代价矩阵，all instance type is DELAYED
		tempLearner.getInitialization(percentage);	
		
		
		// Step 2. iterative active learning
		int[][] recordMatrix = tempLearner.getIterativeLearning();		
		// Step 3. verify the classification result
		int[] correctAndUncorrect = tempLearner.getVerifiedResult();
		
		int[][] returnResult = new int[recordMatrix.length + 2][2];
		returnResult[0] = correctAndUncorrect;
		returnResult[1][0] = tempLearner.totalGuessed;
		returnResult[1][1] = tempLearner.totalTaught;
		for (int i = 2; i < returnResult.length; i ++) {
			returnResult[i] = recordMatrix[i - 2];
		}// Of for i
		
		// record the number of taught instances
		noneTotalTaught = tempLearner.totalTaught;
		
		return returnResult;
	}// Of compareSelf
	
	
	
	/*****
	 * Get predicted class label of target unlabeled instance.
	 * @param paraIndex the index of target unlabeled instance.
	 * @return predictedLabel class label.
	 */
	public int getNoneClassLabel(int paraIndex) {
		// temporary double vector for recording the distances between other labeled instances.
		double[] distanceVector = distanceMatrix[paraIndex];
		
		int count = 0;
		for (int i = 0; i < numInstances(); i ++) {
			if (instanceTypes[i] == TAUGHT) {
				count ++;
			}// Of if
		}// Of for i
		
		// extended a vector to a matrix
		int index = 0;
		double[][] usefulInformation = new double[count][2];
		for (int i = 0; i < numInstances(); i ++) {
			if (instanceTypes[i] == TAUGHT) {
				usefulInformation[index][0] = i;
				usefulInformation[index][1] = distanceVector[i];
				index ++;
			}// Of if
		}// Of for i
		
		// sort
		for (int i = 0; i < usefulInformation.length; i ++) {
			for (int j = i + 1; j < usefulInformation.length; j ++) {
				if (usefulInformation[i][1] > usefulInformation[j][1]) {
					double[] tempArray = usefulInformation[i];
					usefulInformation[i] = usefulInformation[j];
					usefulInformation[j] = tempArray;
				}// Of if
			}// Of for j
		}// Of for i
		
		// find the nearest neighbors
		int[] neighborIndices = new int[kValue];
		int[] vote = new int[numClasses()];
		for (int i = 0; i < neighborIndices.length; i ++) {
			int classLabel = (int)instance(neighborIndices[i]).classValue();
			vote[classLabel] ++;
		}// Of for i
		
		// vote for the predicted class label
		int max = Integer.MIN_VALUE;
		int predictedLabel = -1;
		for (int i = 0; i < vote.length; i ++) {
			if (vote[i] > max) {
				max = vote[i];
				predictedLabel = i;
			}// Of if
		}// Of for i
		
		return predictedLabel;
	}// Of getNoneClassLabel
	
	
	
	/*****
	 * Get kNN classification result.
	 */
	public void getKNNClassification() {
		for (int i = 0; i < numInstances(); i ++) {
			// ignore the instacnes which are taught
			if (instanceTypes[i] == TAUGHT) {
				continue;
			}// Of if
			
			// get the classification result of this unlabeled instance
			int predictedLabel = getNoneClassLabel(i);
			// set value in new added attribute
			instance(i).setValue(0, predictedLabel);
		}// Of for i
	}// Of getKNNClassification
	
	public int[] getNoneVerifiedResult() {
		int tempCorrect = 0;
		int tempUncorrect = 0;
		
		for (int i = 0; i < numInstances(); i ++) {
			int predictedLabel = (int)instance(i).value(0);
			int intrinsicLabel = (int)instance(i).classValue();
			if (predictedLabel == intrinsicLabel) {
				tempCorrect ++;
			} else if (predictedLabel != intrinsicLabel) {
				tempUncorrect ++;
			} else {
				System.out.println("Error occured in getNoneVerifiedResult!");
			}// Of if
		}// Of for i

		int[] tempClassificationResult = new int[2];
		tempClassificationResult[0] = tempCorrect;
		tempClassificationResult[1] = tempUncorrect;
		
		return tempClassificationResult;
	}// Of getNoneVerifiedResult
	
	
	
	/*****
	 * Compare with kNN algorithm.
	 * @param tempLearner the basic decision system.
	 * @return returnResult classification results.
	 */
	public static int[][] compareKNN(Learner tempLearner) {
		// randomly get initial labeled instances
		int labeled = (int)(percentage * tempLearner.numInstances()) + noneTotalTaught;
		double newPercentage = (labeled + 0.0) / tempLearner.numInstances();
		
		// initialization
		tempLearner.getParameterSetting();
		tempLearner.getInitialization(newPercentage);
		
		// get classification results
		tempLearner.getKNNClassification();
		// verify the classification results
		int[] classificationResult = tempLearner.getNoneVerifiedResult();
		int[][] returnResult = new int[1][];
		returnResult[0] = classificationResult;
		
		return returnResult;
	}// Of compareKNN
	
	
	
	/*****
	 * Divide the initial set into training set and testing set.
	 * @param paraPercentage the scale of training set.
	 * @return subsets the set of training set and testing set.
	 * @throws Exception
	 */
	public Learner[] getDividedDataSets(double paraPercentage) throws Exception {
		boolean[] firstInclusionArray = SimpleTool.generateBooleanArrayForDivision(numInstances(), paraPercentage);
		
		
		Learner firstDecisionSystem = new Learner(this);
		firstDecisionSystem.delete(firstInclusionArray);

		
		boolean[] secondInclusionArray = SimpleTool.revertBooleanArray(firstInclusionArray);
		Learner secondDecisionSystem = new Learner(this);
		secondDecisionSystem.delete(secondInclusionArray);

		Learner[] subsets = new Learner[2];
		subsets[0] = firstDecisionSystem;
		subsets[1] = secondDecisionSystem;
		
		return subsets;
	}// Of getDividedDataSets
	
	
	
	/*****
	 * Compare with cost-sensitive decision tree.
	 * @param tempLearner the basic decision system.
	 * @return returnResult the decision-tree-based classification results. 
	 */
	public static int[][] compareDecisionTree(Learner tempLearner) {
		// randomly get initial labeled instances
		int labeled = (int)(percentage * tempLearner.numInstances()) + noneTotalTaught;
		double newPercentage = (labeled + 0.0) / tempLearner.numInstances();
		
		Learner trainingSet = null;
		Learner testingSet = null;
		Learner[] tempSets = null;

		Learner trainingCSId3 = null;
		Learner trainingCSJ48 = null;
		Learner testingCSId3 = null;
		Learner testingCSJ48 = null;
		
		try {
			// divided into training set and testing set
			tempSets = tempLearner.getDividedDataSets(newPercentage);
			trainingSet = tempSets[0];
			testingSet = tempSets[1];
			
			// initial data set
			trainingCSId3 = new Learner(trainingSet);
			trainingCSJ48 = new Learner(trainingSet);
			testingCSId3 = new Learner(testingSet);
			testingCSJ48 = new Learner(testingSet);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error occured in compareDecisionTree!");
		}// Of try
		
		int[] CSId3Result = compareCSId3(trainingCSId3, testingCSId3);
		int[] CSJ48Result = compareCSJ48(trainingCSJ48, testingCSJ48);
		int[][] returnResult = new int[2][];
		returnResult[0] = CSId3Result;
		returnResult[1] = CSJ48Result;
		
		//returnResult[0] = CSId3Result;
		//returnResult[1] = CSJ48Result;
		
		return returnResult;
	}// Of compareDecisionTrees
	
	
	
	/******
	 * Generate Id3-cost-sensitive training set.
	 * @return CSTrainingSet decision system.
	 */
	public Learner getCSId3TrainingSet() {
		Learner CSTrainingSet = null;
		try {
			// Copy the decision system
			CSTrainingSet = new Learner(this);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error occured in getCSId3TrainingSet!");
		}// Of try
		
		// delete instances which missing class label
		CSTrainingSet.deleteWithMissingClass();
		
		//在ID3算法中，选择信息增益最大的属性作为当前的特征对数据集分类。
		Id3 classifier = new Id3();
		
		// set misclassification cost settings
		costMatrixWeka = new CostMatrix(CSTrainingSet.numClasses());
		costMatrixWeka.setCell(0, 0, (0 + 0.0));
		costMatrixWeka.setCell(0, 1, (misCostA + 0.0));
		costMatrixWeka.setCell(1, 0, (misCostB + 0.0));
		costMatrixWeka.setCell(1, 1, (0 + 0.0));
	
//===================================================================		
//		for (int i = 0; i < numClasses(); i ++) {
//			for (int j = 0; j < numClasses(); j ++) {
//				if (i == j) {
//					costMatrixWeka.setCell(i, j, (0 + 0.0));
//				} else {
//					costMatrixWeka.setCell(i, j, (50 + 0.0));
//				}// Of if
//			}// Of for j
//		}// Of for i
//===================================================================
		
		Random random = null;
		if (!(classifier instanceof WeightedInstancesHandler)) {
			random = new Random(1);
		}// Of if
		
		// generate Id3-cost-sensitive set
		try {
			Instances CSDataSet = costMatrixWeka.applyCostMatrix(CSTrainingSet, random);
			CSTrainingSet = new Learner(CSDataSet);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error occured in getCSId3TrainingSet!");
		}// Of try
		
		return CSTrainingSet;
	}// Of getCSId3TrainingSet
	
	/******
	 * Compare classification results with cost sensitive Id3 classifier.
	 * @param trainingSet.
	 * @param testingSet.
	 * @return returnResult.
	 */
	public static int[] compareCSId3(Learner trainingSet, Learner testingSet) {
		// generate Id3-cost-sensitive training set
		Learner CSTrainingSet = trainingSet.getCSId3TrainingSet();
		
		int tempCorrect = 0;
		int tempUncorrect = 0;
		int tempUncorrectA = 0;
		int tempUncorrectB = 0;
		int tempUnclassified = 0;
		
		Id3 classifier = new Id3();
		try {
			// build classifier，buildClassifier方法调用makeTree(Instance data)方法构建分类器，
			//makeTree（Instance data）首先调用computerInfoGain(Instance data,Attribute att)计算每个属性的信息增益
			//选取信息增益最大的属性作为m_Attribute,
			//最后调用Instance[]splitData(Instance data,Attribute att)方法，该方法根据m_Attribute产生与m_Attribute属性取值相当的划分
			//最后在一颗子树上递归的调用makeTree()方法
			classifier.buildClassifier(CSTrainingSet);
			//测试集
			for (int i = 0; i < testingSet.numInstances(); i ++) {
				// get predicted class label
				int predictedLabel = (int)classifier.classifyInstance(testingSet.instance(i));
				// get actual class label
				
				int trueLabel = (int)testingSet.instance(i).classValue();		
				//tempUnclassified==0???
				if (trueLabel < 0 || trueLabel >= testingSet.numClasses()) {
					tempUnclassified ++;
				} else {
					if (trueLabel == predictedLabel) {
						tempCorrect ++;
					} else {
						if (trueLabel == 0) {
							tempUncorrectA ++;
						} else if (trueLabel == 1) {
							tempUncorrectB ++;
						} else {
							System.out.println("Error occured in compareCSId3!");
						}// Of if
						
						tempUncorrect ++;
					}// Of if
				}// Of if
			}// Of for iError occured in getVerifiedResult!
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error occured in compareCSId3!");
		}// Of try
		
		int[] returnResult = new int[5];
		returnResult[0] = tempCorrect;
		returnResult[1] = tempUncorrect;
		returnResult[2] = tempUncorrectA;
		returnResult[3] = tempUncorrectB;
		returnResult[4] = tempUnclassified;
		
		
		return returnResult;
	}// Of compareCSId3
	
	
	
	/*****
	 * Generate J48-cost-sensitive training set.
	 * @return CSTrainingSet decision system.
	 */
	public Learner getCSJ48TrainingSet() {
		Learner CSTrainingSet = null;
		try {
			// Copy the decision system
			CSTrainingSet = new Learner(this);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error occured in getCSJ48TrainingSet!");
		}// Of try
		
		// delete instances which missing class label
		CSTrainingSet.deleteWithMissingClass();
		J48 classifier = new J48();
		
		// set misclassification cost settings
		costMatrixWeka = new CostMatrix(CSTrainingSet.numClasses());
		costMatrixWeka.setCell(0, 0, 0 + 0.0);
		costMatrixWeka.setCell(0, 1, misCostA + 0.0);
		costMatrixWeka.setCell(1, 0, misCostB + 0.0);
		costMatrixWeka.setCell(1, 1, 0 + 0.0);
		
//===================================================================		
//		for (int i = 0; i < numClasses(); i ++) {
//			for (int j = 0; j < numClasses(); j ++) {
//				if (i == j) {
//					costMatrixWeka.setCell(i, j, 0);
//				} else {
//					costMatrixWeka.setCell(i, j, 50);
//				}// Of if
//			}// Of for j
//		}// Of for i
//===================================================================
		
		Random random = null;
		if (!(classifier instanceof WeightedInstancesHandler)) {
			random = new Random(1);
		}// Of if
		
		try {
			// generate J48-cost-sensitive set
			Instances CSDataSet = costMatrixWeka.applyCostMatrix(CSTrainingSet, random);
			CSTrainingSet = new Learner(CSDataSet);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error occured in getCSId3TrainingSet!");
		}// Of try
		
		return CSTrainingSet;
	}// Of getCSJ48TrainingSet
	
	public static int[] compareCSJ48(Learner trainingSet, Learner testingSet) {
		int tempCorrect = 0;
		int tempUncorrect = 0;
		int tempUncorrectA = 0;
		int tempUncorrectB = 0;
		int tempUnclassified = 0;
		
		J48 classifier = new J48();
		try {
			// generate J48-cost-sensitive training set
			classifier.buildClassifier(trainingSet);
		
			for (int i = 0; i < testingSet.numInstances(); i ++) {
				// get predicted class label
				int predictedLabel = (int)classifier.classifyInstance(testingSet.instance(i));
				// get actual class label
				int trueLabel = (int)testingSet.instance(i).classValue();
				if (trueLabel < 0 || trueLabel >= testingSet.numClasses()) {
					tempUnclassified ++;
				} else {
					if (trueLabel == predictedLabel) {
						tempCorrect ++;
					} else {
						if (trueLabel == 0) {
							tempUncorrectA ++;
						} else if (trueLabel == 1) {
							tempUncorrectB ++;
						} else {
							System.out.println("Error occured in compareCSJ48!");
						}// Of if
						
						tempUncorrect ++;
					}// Of if
				}// Of if
			}// Of for i
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error occured in compareCSJ48!");
		}// Of try
		
		int[] returnResult = new int[5];
		returnResult[0] = tempCorrect;
		returnResult[1] = tempUncorrect;
		returnResult[2] = tempUncorrectA;
		returnResult[3] = tempUncorrectB;
		returnResult[4] = tempUnclassified;
		
		return returnResult;
	}// Of compareCSJ48
}// Of class Learner