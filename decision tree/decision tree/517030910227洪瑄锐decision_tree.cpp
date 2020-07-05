/*
517030910227 洪瑄锐

该决策树应用的算法为ID3和C4.5
在函数build_decision_tree中选择两个语句中的一个注释掉即可更换算法
该决策树应用的数据集为判断乳腺癌为良性（数据集中以2代表）或为恶性（数据集中以4为代表）
该数据集的属性共有10个，训练和测试时选择了属性4-10，每个属性的值为1-10
运行代码输出决策树，输出将训练集放入决策树的成功率，输出将测试集放入决策树的成功率
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <sstream>

using namespace std;

vector<vector<unsigned>> trains;//训练集内容
vector<vector<unsigned>> tests;

//属性名称，作用于打印决策树
string attribute_names[] = { "Uniformity of Cell Shape", "Marginal Adhesion", "Single Epithelial Cell Size","Bare Nuclei","Bland Chromatin","Normal Nucleoli","Mitoses" };

//数据集中各个属性的所有可能取值，与attribute_names[]中的元素按顺序对应
unsigned attribute_values[] = { 1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,2,4 };

//将所有属性和分类的可能取值与0-71一 一对应，读取数据集时将数据转化为对应的attribute_number
unsigned attribute_number[] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71 };


//将数据读入训练集和测试集，并转换为对应的数字，存入动态数组trains和tests

void read_file()
{
	vector<unsigned> s;//暂存数据集中的每一行

	unsigned shape;//10种属性和一种类别，其中3种属性忽略
	unsigned adhesion;
	unsigned size;
	unsigned nuclei;
	unsigned chromatin;
	unsigned nucleoli;
	unsigned mitoses;
	unsigned classes;
	unsigned codeNumber;
	unsigned thickness;
	unsigned cellSize;

	//训练集
	fstream is;
	is.open("breast-cancer.data");

	if (is)
	{
		while (!is.eof())
		{
			is >> codeNumber >> thickness >> cellSize >> shape >> adhesion >> size >> nuclei >> chromatin >> nucleoli >> mitoses >> classes;
			//转换为对应数字
			s.push_back(shape - 1);
			s.push_back(adhesion + 9);
			s.push_back(size + 19);
			s.push_back(nuclei + 29);
			s.push_back(chromatin + 39);
			s.push_back(nucleoli + 49);
			s.push_back(mitoses + 59);
			if (classes == 2)
				s.push_back(70);
			else
				s.push_back(71);

			trains.push_back(s);
			s.clear();
		}
	}
	is.close();

	//测试集
	fstream tis;
	tis.open("test.data");
	if (tis)
	{
		while (!tis.eof())
		{
			tis >> codeNumber >> thickness >> cellSize >> shape >> adhesion >> size >> nuclei >> chromatin >> nucleoli >> mitoses >> classes;

			s.push_back(shape - 1);
			s.push_back(adhesion + 9);
			s.push_back(size + 19);
			s.push_back(nuclei + 29);
			s.push_back(chromatin + 39);
			s.push_back(nucleoli + 49);
			s.push_back(mitoses + 59);
			if (classes == 2)
				s.push_back(70);
			else
				s.push_back(71);
			tests.push_back(s);
			s.clear();
		}
	}
	tis.close();
}


//计算一个数值得以2为底的对数

double log2(double n)
{
	return log10(n) / log10(2.0);
}

//将vector中重复元素合并，只保留一个

template <typename T>
vector<T> unique(vector<T> vals)
{
	vector<T> unique_vals;
	vector<T>::iterator itr;
	vector<T>::iterator subitr;

	int flag = 0;
	while (!vals.empty())
	{
		unique_vals.push_back(vals[0]);
		itr = vals.begin();
		subitr = unique_vals.begin() + flag;
		while (itr != vals.end())
		{
			if (*subitr == *itr)
				itr = vals.erase(itr);
			else
				itr++;
		}
		flag++;
	}
	return unique_vals;
}

//根据属性的取值，计算该属性的熵

double compute_entropy(vector<unsigned> v) //输入为不同组数据中某种属性或分类的取值集合，例如classes：（2，2，4，4，2，2，4，4，2......）
{
	vector<unsigned> unique_v;
	unique_v = unique(v);//去掉重复元素，在本数据集中即（2，4）

	vector<unsigned>::iterator itr;
	itr = unique_v.begin();

	double entropy = 0.0;
	auto total = v.size();
	while (itr != unique_v.end())//遍历classes的每种取值，计算 概率*log2（概率），取负后加入entropy
	{
		double cnt = count(v.begin(), v.end(), *itr);//计算每种classes取值在集合中的个数
		entropy -= cnt / total * log2(cnt / total);
		itr++;
	}
	return entropy;
}

//计算数据集中所有属性的信息增益
vector<double> compute_gain(vector<vector<unsigned> > truths)//输入为训练集
{
	vector<double> gain(truths[0].size() - 1, 0);//存储各个属性的信息增益，规模数为truths的size-1是因为训练集各个数据的最后一列为classes，无需计算
	vector<unsigned> attribute_vals;//用来针对每个属性存储其取值，得到该属性在各个数据中的取值集合
	vector<unsigned> labels;//用来存储classes的取值
	for (unsigned j = 0; j < truths.size(); j++)//将classes取值放入labels
	{
		labels.push_back(truths[j].back());
	}

	for (unsigned i = 0; i < truths[0].size() - 1; i++)//遍历每一种属性
	{
		for (unsigned j = 0; j < truths.size(); j++)//将每个数据j中的第i种属性的取值放入attribute_vals
			attribute_vals.push_back(truths[j][i]);

		vector<unsigned> unique_vals = unique(attribute_vals);//去重，得到第i种属性的所有可能取值
		vector<unsigned>::iterator itr = unique_vals.begin();
		vector<unsigned> subset;//subset存储对应于数据的第i种属性为某个取值时的classes集合
		while (itr != unique_vals.end())//遍历第i种属性的所有可能取值
		{
			for (unsigned k = 0; k < truths.size(); k++)
			{
				if (*itr == attribute_vals[k])
				{
					subset.push_back(truths[k].back());//将数据truth[k]的第i种属性取*itr时的classes放入subset
				}
			}
			double A = (double)subset.size();
			gain[i] += A / truths.size() * compute_entropy(subset);//得到第i种属性取*itr的概率以及该属性取值对应classes的信息熵，加到H(D|A)中
			itr++;//计算第i种属性的下一个取值
			subset.clear();
		}
		gain[i] = compute_entropy(labels) - gain[i];//g(D,A)=H(D)-H(D|A)
		attribute_vals.clear();
	}
	return gain;
}

//计算数据集中所有属性的信息增益率

vector<double> compute_gain_ratio(vector<vector<unsigned> > truths)
{
	vector<double> gain = compute_gain(truths);//获得每种属性的信息增益熵
	vector<double> entropies;
	vector<double> gain_ratio;

	for (unsigned i = 0; i < truths[0].size() - 1; i++)//最后一列是类别标签，不需要计算信息增益率
	{
		vector<unsigned> attribute_vals(truths.size(), 0);
		for (unsigned j = 0; j < truths.size(); j++)
		{
			attribute_vals[j] = truths[j][i];
		}
		double current_entropy = compute_entropy(attribute_vals);
		if (current_entropy)
		{
			gain_ratio.push_back(gain[i] / current_entropy);
		}
		else
			gain_ratio.push_back(0.0);

	}
	return gain_ratio;
}

//找出数据集中最多的类别标签

template <typename T>
T find_most_common_label(vector<vector<T> > data)
{
	vector<T> labels;
	for (unsigned i = 0; i < data.size(); i++)
	{
		labels.push_back(data[i].back());
	}
	vector<T>::iterator itr = labels.begin();
	T most_common_label;
	unsigned most_counter = 0;
	while (itr != labels.end())
	{
		unsigned current_counter = count(labels.begin(), labels.end(), *itr);
		if (current_counter > most_counter)
		{
			most_common_label = *itr;
			most_counter = current_counter;
		}
		itr++;
	}
	return most_common_label;
}


//根据属性，找出该属性可能的取值

template <typename T>
vector<T> find_attribute_values(T attribute, vector<vector<T> > data)
{
	vector<T> values;
	for (unsigned i = 0; i < data.size(); i++)
	{
		values.push_back(data[i][attribute]);
	}
	return unique(values);
}

/*
在构建决策树的过程中，如果某一属性已经考察过了
那么就从数据集中去掉这一属性，此处不是真正意义
上的去掉，而是将考虑过的属性全部标记为100
*/

template <typename T>
vector<vector<T> > drop_one_attribute(T attribute, vector<vector<T> > data)
{
	vector<vector<T> > new_data(data.size(), vector<T>(data[0].size() - 1, 0));
	for (unsigned i = 0; i < data.size(); i++)
	{
		data[i][attribute] = 100;
	}
	return data;
}


struct Tree {
	unsigned root;//节点属性值
	vector<unsigned> branches;//节点可能取值
	vector<Tree> children; //孩子节点
};

//递归构建决策树
void build_decision_tree(vector<vector<unsigned> > examples, Tree &tree)
{
	//递归终止条件：判断所有实例是否都属于同一类，如果是，则决策树是单节点
	vector<unsigned> labels(examples.size(), 0);
	for (unsigned i = 0; i < examples.size(); i++)
	{
		labels[i] = examples[i].back();
	}
	if (unique(labels).size() == 1)
	{
		tree.root = labels[0];
		return;
	}

	//递归终止条件：判断是否还有剩余的属性没有考虑，如果所有属性都已经考虑过了，
	//那么此时属性数量为0，将训练集中最多的类别标记作为该节点的类别标记
	if (count(examples[0].begin(), examples[0].end(), 110) == examples[0].size() - 1)//只剩下一列类别标记
	{
		tree.root = find_most_common_label(examples);
		return;
	}
	//计算信息增益，选择信息增益最大的属性作为根节点,并找出该节点的所有取值
	
	//vector<double> standard = compute_gain(examples);

	//要是采用C4.5，将上面一行注释掉，把下面一行的注释去掉即可
	vector<double> standard = compute_gain_ratio(examples);

	tree.root = 0;
	for (unsigned i = 0; i < standard.size(); i++)
	{
		if (standard[i] >= standard[tree.root] && examples[0][i] != 100)
			tree.root = i;
	}

	tree.branches = find_attribute_values(tree.root, examples);
	//根据节点的取值，将examples分成若干子集
	vector<vector<unsigned> > new_examples = drop_one_attribute(tree.root, examples);
	vector<vector<unsigned> > subset;
	for (unsigned i = 0; i < tree.branches.size(); i++)
	{
		for (unsigned j = 0; j < examples.size(); j++)
		{
			for (unsigned k = 0; k < examples[0].size(); k++)
			{
				if (tree.branches[i] == examples[j][k])
					subset.push_back(new_examples[j]);
			}
		}
	    //对每一个子集递归调用build_decision_tree()函数
		Tree new_tree;
		build_decision_tree(subset, new_tree);
		tree.children.push_back(new_tree);
		subset.clear();
	}
}

//打印决策树

void print_decision_tree(Tree tree, unsigned depth)
{
	for (unsigned d = 0; d < depth; d++) cout << "\t";
	if (!tree.branches.empty()) //不是叶子节点
	{
		cout << attribute_names[tree.root] << endl;

		for (unsigned i = 0; i < tree.branches.size(); i++)
		{
			for (unsigned d = 0; d < depth + 1; d++) cout << "\t";
			cout << attribute_values[tree.branches[i]] << endl;
			print_decision_tree(tree.children[i], depth + 2);
		}
	}
	else //是叶子节点
	{
		cout << attribute_values[tree.root] << endl;
	}

}

//根据数据属性取值预测类别

unsigned classify_tree(Tree tree, vector<unsigned> test)
{
	if (tree.branches.empty())//是叶子节点
	{
		return tree.root;
	}
	else
	{
		for (unsigned i = 0; i < tree.branches.size(); i++)
		{
			if (test[tree.root] == tree.branches[i])
			{
				return classify_tree(tree.children[i], test);
			}
		}
	}
}

//检测决策树学习效果和推广性能

void test_tree(Tree tree)
{
	unsigned num_right_test = 0;
	unsigned num_right_train = 0;

	unsigned result;

	double right_ratio_test;
	double right_ratio_train;

	//训练集，得出学习效果
	for (unsigned i = 0; i < trains.size(); i++)
	{
		result = classify_tree(tree, trains[i]);
		if (result == trains[i].back())
			num_right_train++;
	}
	right_ratio_train = double(num_right_train) / trains.size();
	cout << endl;
	cout << "训练集成功率为" << right_ratio_train;

	//测试集，得出推广性能
	for (unsigned i = 0; i < tests.size(); i++)
	{
		result = classify_tree(tree, tests[i]);
		if (result == tests[i].back())
			num_right_test++;
	}
	right_ratio_test = double(num_right_test) / tests.size();
	cout << endl;
	cout << "测试集成功率为" << right_ratio_test;
}



int main()
{

	read_file();//读取训练集的数据入动态数组trains,读取测试集的数据入动态数组tests

	Tree tree;//创建一棵决策树
	build_decision_tree(trains, tree);

	print_decision_tree(tree, 0);//打印决策树

	test_tree(tree);//将训练集数据和测试集数据分别放入决策树进行测试，输出各自的预测成功率，得到学习效果和推广性能

	return 0;
}