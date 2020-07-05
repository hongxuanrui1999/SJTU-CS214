/*
517030910227 ��u��

�þ�����Ӧ�õ��㷨ΪID3��C4.5
�ں���build_decision_tree��ѡ����������е�һ��ע�͵����ɸ����㷨
�þ�����Ӧ�õ����ݼ�Ϊ�ж����ٰ�Ϊ���ԣ����ݼ�����2������Ϊ���ԣ����ݼ�����4Ϊ����
�����ݼ������Թ���10����ѵ���Ͳ���ʱѡ��������4-10��ÿ�����Ե�ֵΪ1-10
���д�������������������ѵ��������������ĳɹ��ʣ���������Լ�����������ĳɹ���
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

vector<vector<unsigned>> trains;//ѵ��������
vector<vector<unsigned>> tests;

//�������ƣ������ڴ�ӡ������
string attribute_names[] = { "Uniformity of Cell Shape", "Marginal Adhesion", "Single Epithelial Cell Size","Bare Nuclei","Bland Chromatin","Normal Nucleoli","Mitoses" };

//���ݼ��и������Ե����п���ȡֵ����attribute_names[]�е�Ԫ�ذ�˳���Ӧ
unsigned attribute_values[] = { 1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,2,4 };

//���������Ժͷ���Ŀ���ȡֵ��0-71һ һ��Ӧ����ȡ���ݼ�ʱ������ת��Ϊ��Ӧ��attribute_number
unsigned attribute_number[] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71 };


//�����ݶ���ѵ�����Ͳ��Լ�����ת��Ϊ��Ӧ�����֣����붯̬����trains��tests

void read_file()
{
	vector<unsigned> s;//�ݴ����ݼ��е�ÿһ��

	unsigned shape;//10�����Ժ�һ���������3�����Ժ���
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

	//ѵ����
	fstream is;
	is.open("breast-cancer.data");

	if (is)
	{
		while (!is.eof())
		{
			is >> codeNumber >> thickness >> cellSize >> shape >> adhesion >> size >> nuclei >> chromatin >> nucleoli >> mitoses >> classes;
			//ת��Ϊ��Ӧ����
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

	//���Լ�
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


//����һ����ֵ����2Ϊ�׵Ķ���

double log2(double n)
{
	return log10(n) / log10(2.0);
}

//��vector���ظ�Ԫ�غϲ���ֻ����һ��

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

//�������Ե�ȡֵ����������Ե���

double compute_entropy(vector<unsigned> v) //����Ϊ��ͬ��������ĳ�����Ի�����ȡֵ���ϣ�����classes����2��2��4��4��2��2��4��4��2......��
{
	vector<unsigned> unique_v;
	unique_v = unique(v);//ȥ���ظ�Ԫ�أ��ڱ����ݼ��м���2��4��

	vector<unsigned>::iterator itr;
	itr = unique_v.begin();

	double entropy = 0.0;
	auto total = v.size();
	while (itr != unique_v.end())//����classes��ÿ��ȡֵ������ ����*log2�����ʣ���ȡ�������entropy
	{
		double cnt = count(v.begin(), v.end(), *itr);//����ÿ��classesȡֵ�ڼ����еĸ���
		entropy -= cnt / total * log2(cnt / total);
		itr++;
	}
	return entropy;
}

//�������ݼ����������Ե���Ϣ����
vector<double> compute_gain(vector<vector<unsigned> > truths)//����Ϊѵ����
{
	vector<double> gain(truths[0].size() - 1, 0);//�洢�������Ե���Ϣ���棬��ģ��Ϊtruths��size-1����Ϊѵ�����������ݵ����һ��Ϊclasses���������
	vector<unsigned> attribute_vals;//�������ÿ�����Դ洢��ȡֵ���õ��������ڸ��������е�ȡֵ����
	vector<unsigned> labels;//�����洢classes��ȡֵ
	for (unsigned j = 0; j < truths.size(); j++)//��classesȡֵ����labels
	{
		labels.push_back(truths[j].back());
	}

	for (unsigned i = 0; i < truths[0].size() - 1; i++)//����ÿһ������
	{
		for (unsigned j = 0; j < truths.size(); j++)//��ÿ������j�еĵ�i�����Ե�ȡֵ����attribute_vals
			attribute_vals.push_back(truths[j][i]);

		vector<unsigned> unique_vals = unique(attribute_vals);//ȥ�أ��õ���i�����Ե����п���ȡֵ
		vector<unsigned>::iterator itr = unique_vals.begin();
		vector<unsigned> subset;//subset�洢��Ӧ�����ݵĵ�i������Ϊĳ��ȡֵʱ��classes����
		while (itr != unique_vals.end())//������i�����Ե����п���ȡֵ
		{
			for (unsigned k = 0; k < truths.size(); k++)
			{
				if (*itr == attribute_vals[k])
				{
					subset.push_back(truths[k].back());//������truth[k]�ĵ�i������ȡ*itrʱ��classes����subset
				}
			}
			double A = (double)subset.size();
			gain[i] += A / truths.size() * compute_entropy(subset);//�õ���i������ȡ*itr�ĸ����Լ�������ȡֵ��Ӧclasses����Ϣ�أ��ӵ�H(D|A)��
			itr++;//�����i�����Ե���һ��ȡֵ
			subset.clear();
		}
		gain[i] = compute_entropy(labels) - gain[i];//g(D,A)=H(D)-H(D|A)
		attribute_vals.clear();
	}
	return gain;
}

//�������ݼ����������Ե���Ϣ������

vector<double> compute_gain_ratio(vector<vector<unsigned> > truths)
{
	vector<double> gain = compute_gain(truths);//���ÿ�����Ե���Ϣ������
	vector<double> entropies;
	vector<double> gain_ratio;

	for (unsigned i = 0; i < truths[0].size() - 1; i++)//���һ��������ǩ������Ҫ������Ϣ������
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

//�ҳ����ݼ�����������ǩ

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


//�������ԣ��ҳ������Կ��ܵ�ȡֵ

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
�ڹ����������Ĺ����У����ĳһ�����Ѿ��������
��ô�ʹ����ݼ���ȥ����һ���ԣ��˴�������������
�ϵ�ȥ�������ǽ����ǹ�������ȫ�����Ϊ100
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
	unsigned root;//�ڵ�����ֵ
	vector<unsigned> branches;//�ڵ����ȡֵ
	vector<Tree> children; //���ӽڵ�
};

//�ݹ鹹��������
void build_decision_tree(vector<vector<unsigned> > examples, Tree &tree)
{
	//�ݹ���ֹ�������ж�����ʵ���Ƿ�����ͬһ�࣬����ǣ���������ǵ��ڵ�
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

	//�ݹ���ֹ�������ж��Ƿ���ʣ�������û�п��ǣ�����������Զ��Ѿ����ǹ��ˣ�
	//��ô��ʱ��������Ϊ0����ѵ�����������������Ϊ�ýڵ�������
	if (count(examples[0].begin(), examples[0].end(), 110) == examples[0].size() - 1)//ֻʣ��һ�������
	{
		tree.root = find_most_common_label(examples);
		return;
	}
	//������Ϣ���棬ѡ����Ϣ��������������Ϊ���ڵ�,���ҳ��ýڵ������ȡֵ
	
	//vector<double> standard = compute_gain(examples);

	//Ҫ�ǲ���C4.5��������һ��ע�͵���������һ�е�ע��ȥ������
	vector<double> standard = compute_gain_ratio(examples);

	tree.root = 0;
	for (unsigned i = 0; i < standard.size(); i++)
	{
		if (standard[i] >= standard[tree.root] && examples[0][i] != 100)
			tree.root = i;
	}

	tree.branches = find_attribute_values(tree.root, examples);
	//���ݽڵ��ȡֵ����examples�ֳ������Ӽ�
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
	    //��ÿһ���Ӽ��ݹ����build_decision_tree()����
		Tree new_tree;
		build_decision_tree(subset, new_tree);
		tree.children.push_back(new_tree);
		subset.clear();
	}
}

//��ӡ������

void print_decision_tree(Tree tree, unsigned depth)
{
	for (unsigned d = 0; d < depth; d++) cout << "\t";
	if (!tree.branches.empty()) //����Ҷ�ӽڵ�
	{
		cout << attribute_names[tree.root] << endl;

		for (unsigned i = 0; i < tree.branches.size(); i++)
		{
			for (unsigned d = 0; d < depth + 1; d++) cout << "\t";
			cout << attribute_values[tree.branches[i]] << endl;
			print_decision_tree(tree.children[i], depth + 2);
		}
	}
	else //��Ҷ�ӽڵ�
	{
		cout << attribute_values[tree.root] << endl;
	}

}

//������������ȡֵԤ�����

unsigned classify_tree(Tree tree, vector<unsigned> test)
{
	if (tree.branches.empty())//��Ҷ�ӽڵ�
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

//��������ѧϰЧ�����ƹ�����

void test_tree(Tree tree)
{
	unsigned num_right_test = 0;
	unsigned num_right_train = 0;

	unsigned result;

	double right_ratio_test;
	double right_ratio_train;

	//ѵ�������ó�ѧϰЧ��
	for (unsigned i = 0; i < trains.size(); i++)
	{
		result = classify_tree(tree, trains[i]);
		if (result == trains[i].back())
			num_right_train++;
	}
	right_ratio_train = double(num_right_train) / trains.size();
	cout << endl;
	cout << "ѵ�����ɹ���Ϊ" << right_ratio_train;

	//���Լ����ó��ƹ�����
	for (unsigned i = 0; i < tests.size(); i++)
	{
		result = classify_tree(tree, tests[i]);
		if (result == tests[i].back())
			num_right_test++;
	}
	right_ratio_test = double(num_right_test) / tests.size();
	cout << endl;
	cout << "���Լ��ɹ���Ϊ" << right_ratio_test;
}



int main()
{

	read_file();//��ȡѵ�����������붯̬����trains,��ȡ���Լ��������붯̬����tests

	Tree tree;//����һ�þ�����
	build_decision_tree(trains, tree);

	print_decision_tree(tree, 0);//��ӡ������

	test_tree(tree);//��ѵ�������ݺͲ��Լ����ݷֱ������������в��ԣ�������Ե�Ԥ��ɹ��ʣ��õ�ѧϰЧ�����ƹ�����

	return 0;
}