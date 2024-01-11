#include <iostream>
#include <cmath>
#include <random>
#include <cstring>
#include <cstdlib>
#include <vector>


using namespace std;

class SomeObject {
public:
    double num;
    char line;

    SomeObject() {
        num = 1.0;
        line = '\0';
    }

    SomeObject(double num, char line) : num(num), line(line) {}

    friend ostream& operator<<(ostream& os, const SomeObject& so) {
        os << "Num: " << so.num << ", Line: " << so.line;
        return os;
    }
};

template <typename T>
int compare(T a, T b)
{
    if (a < b) {
        return -1;
    }
    else if (a > b) {
        return 1;
    }
    else {
        return 0;
    }

};

int numCompare(SomeObject* a, SomeObject* b) {
    return compare(a->num, b->num);
}

double objectKey(SomeObject* so) {
    return so->num;
}

SomeObject randomObject() {

    static random_device rd;
    static default_random_engine gen(1/*rd()*/);

    uniform_real_distribution<double> numDistribution;
    double randomNum = numDistribution(gen);

   /* uniform_int_distribution<int> lineLengthDistribution(5, 15);
    int lineLength = lineLengthDistribution(gen);*/

    uniform_int_distribution<int> charDistribution('a', 'z');
    string randomLine;
    //for (int i = 0; i < lineLength; ++i) {
        randomLine += charDistribution(gen);
    //}

    return { randomNum, randomLine.front() };
}

template<typename T>
class MaxHeap {
public:
    T* array;
    int size;
    int capacity;

    MaxHeap() : size(0), capacity(1) {
        array = new T[capacity];
    }

    MaxHeap(T* array, int capacity, int (*compar)(T, T), bool flag) {
        this->array = array;
        this->size = capacity;
        this->capacity = capacity;

        if (flag) {
            for (int i = 1; i < size; i++) {
                heapifyUp(i, compar);
            }
        }
        else {
            for (int i = (size / 2) - 1; i >= 0; i--) {
                heapifyDown(i, compar);
            }
        }
    }

    ~MaxHeap() {
        delete[] array;
        array = nullptr;
        size = 0;
        capacity = 0;
    }

    int parent(T number) { return (number - 1) / 2; }

    int left(T number) { return (2 * number + 1); }

    int right(T number) { return (2 * number + 2); }

    void swap(T& x, T& y) {
        T temp = x;
        x = y;
        y = temp;
    }

    void heapifyUp(int index, int(*compare)(T, T)) {
        int p = (index - 1) / 2;//parent(index);
        if (index != 0 && compare(array[index], array[p]) > 0) {
            swap(array[index], array[p]);
            heapifyUp(p, compare);
        }
    }

    void heapifyDown(int index, int(*compare)(T, T)) {
        int l = 2 * index + 1;//left(index);
        int r = 2 * index + 2;//right(index);
        int max = index;

        if (r < size && compare(array[r], array[max]) > 0) {
            max = r;
        }
        if (l < size && compare(array[l], array[max]) > 0) {
            max = l;
        }
        if (max != index) {
            swap(array[index], array[max]);
            heapifyDown(max, compare);
        }
    }

    void sort(int(*compare)(T, T)) {
        for (int i = size - 1; i >= 0; i--) {
            swap(array[0], array[i]);
            size--;
            heapifyDown(0, compare);
        }
        size = capacity;
    }

    void deleteHeap() {
        delete[] array;
        array = nullptr;
        size = 0;
        capacity = 0;
    }
};

void countingSortInteger(int* array, int size, int m) {
    int* count = new int[m]();
    int* result = new int[size];

    for (int i = 0; i < size; i++) {
        count[array[i]]++;
    }

    for (int i = 1; i < m; i++) {
        count[i] += count[i - 1];
    }

    for (int i = size - 1; i >= 0; i--) {
        result[count[array[i]] - 1] = array[i];
        count[array[i]]--;
    }

    for (int i = 0; i < size; i++) {
        array[i] = result[i];
    }

    delete[] count;
    delete[] result;
}

void bubbleSort(int* array, int size) {
    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - i - 1; j++) {
            if (array[j] > array[j + 1]) {
                int temp = array[j];
                array[j] = array[j + 1];
                array[j + 1] = temp;
            }
        }
    }
}

void bucketSortInt(int* array, int n, int m) {
    if (array == nullptr || n <= 1) {
        return;
    }

    int maxVal = array[0];
    int minVal = array[0];

    for (int i = 1; i < n; ++i) {
        if (array[i] > maxVal) {
            maxVal = array[i];
        }
        if (array[i] < minVal) {
            minVal = array[i];
        }
    }

    int numBuckets = m;
    vector<vector<int>> buckets(numBuckets);

    for (int i = 0; i < n; ++i) {
        int index = static_cast<int>((array[i] - minVal) * numBuckets / (maxVal - minVal + 1));
        buckets[index].push_back(array[i]);
    }

    for (int i = 0; i < numBuckets; ++i) {
        int bucketSize = static_cast<int>(buckets[i].size());
        if (bucketSize > 1) {
            bubbleSort(buckets[i].data(), bucketSize);
        }
    }

    int index = 0;
    for (int i = 0; i < numBuckets; ++i) {
        for (int j = 0; j < static_cast<int>(buckets[i].size()); ++j) {
            array[index++] = buckets[i][j];
        }
    }
}

template <typename T>
using KeyFunction = double (*)(T);

template <typename T>
using ComparatorFunction = int (*)(T, T);

template <typename T>
void insertionSort(vector<T>& array, ComparatorFunction<T> comparator) {
    for (int i = 1; i < array.size(); ++i) {
        T key = array[i];
        int j = i - 1;
        while (j >= 0 && comparator(array[j], key) > 0) {
            array[j + 1] = array[j];
            --j;
        }
        array[j + 1] = key;
    }
}

template <typename T>
void bucketSortObject(T* array, int n, double m, KeyFunction<T> keyFunction, ComparatorFunction<T> comparator) {
    int num_buckets = n;

    vector<vector<T>> buckets(num_buckets);


    for (int i = 0; i < n; i++) {
        int index = static_cast<int>(keyFunction(array[i]) * n);
        buckets[index].push_back(array[i]);
    }

    for (int i = 0; i < num_buckets; i++) {
        int bucketSize = static_cast<int>(buckets[i].size());
        if (buckets[i].size() > 1) {
            insertionSort(buckets[i], comparator);
        }
    }

    int idx = 0;
    for (int i = 0; i < num_buckets; i++) {
        for (int j = 0; j < static_cast<int>(buckets[i].size()); ++j) {
              array[idx++] = buckets[i][j];
        }
    }
}

template<typename T>
void printArrayInteger(T* array, int size) {
    for (int i = 0; i < size; i++) {
        cout << array[i] << " ";
    }
    cout << endl;
}

void printArrayObjects(SomeObject** so, int n) {
    for (int i = 0; i < n; i++) {
        cout << *so[i] << endl;
    }
}

int randomValue() {
    random_device rd;
    default_random_engine dfe(rd());
    uniform_int_distribution<int> number(1, 200);
    return number(dfe);
}

//int main() {
//    const int MAX_ORDER = 6;
//    const int m = (int) pow (10 , 7);
//
//    for ( int o = 1; o <= MAX_ORDER ; o ++) {
//
//        const int n = (int) pow(10, o);
//        int *array1 = new int[n];
//        for (int i = 0; i < n; i++) {
//            int rand_val = randomValue();
//            array1[i] = rand_val;
//        }
//        //cout << "Main array: ";
//        //printArrayInteger(array1, n);
//
//        int* array2 = new int [n];
//        int* array3 = new int [n];
//        memcpy(array2, array1, n*sizeof(int));
//        memcpy(array3, array1, n*sizeof(int));
//
//        clock_t t1 = clock () ;
//        countingSortInteger(array1, n, m);
//        clock_t t2 = clock () ;
//        double time = static_cast<double>(t2 - t1) / CLOCKS_PER_SEC;
//        cout << "Time to sort by counting sort " << n << " element: " << time << " seconds" << endl;
//        //printArrayInteger(array1, n);
//
//
//        t1 = clock ();
//        auto* bh = new MaxHeap<int>(array2, n, compare, false);
//        bh->sort(compare);
//        t2 = clock();
//        time = static_cast<double>(t2 - t1) / CLOCKS_PER_SEC;
//        cout << "Time to sort by heap sort " << n << " element: " << time << " seconds" << endl;
//        //printArrayInteger(array2, n);
//
//        t1 = clock ();
//        bucketSortInt(array3, n, m);
//        t2 = clock();
//        time = static_cast<double>(t2 - t1) / CLOCKS_PER_SEC;
//        cout << "Time to sort by bucket sort " << n << " element: " << time << " seconds" << endl;
//        //printArrayInteger(array3, n);
//        cout << endl;
//
//        delete[] array1;
//        delete[] array2;
//        delete[] array3;
//    }
//
//    return 0;
//}

int main() {
    const int MAX_ORDER = 7;

    for (int o = 1; o <= MAX_ORDER; o++) {
        const int n = (int)pow(10, o);
        SomeObject** array1 = new SomeObject * [n];
        for (int i = 0; i < n; i++) {
            array1[i] = new SomeObject(randomObject());
        }

        //cout << "Main array:" << endl;
        //printArrayObjects(array1, n);
        //cout << endl;

        SomeObject** array2 = new SomeObject * [n];
        memcpy(array2, array1, n * sizeof(SomeObject*));

        clock_t t1 = clock();
        MaxHeap<SomeObject*>* bh = new MaxHeap<SomeObject*>(array1, n, numCompare, false);
        bh->sort(numCompare);
        clock_t t2 = clock();
        double time = static_cast<double>(t2 - t1) / CLOCKS_PER_SEC;
        cout << "Time to sort by heap sort " << n << " objects: " << time << " seconds" << endl;
        //printArrayObjects(array1, n);

        t1 = clock();
        bucketSortObject<SomeObject*>(array2, n, 1.0, objectKey, numCompare);
        t2 = clock();
        time = static_cast<double>(t2 - t1) / CLOCKS_PER_SEC;
        cout << "Time to sort by bucket sort " << n << " objects: " << time << " seconds" << endl;
        //printArrayObjects(array2, n);

        int trueCounter = 0;
        for (int i = 0; i < n; i++) {
            if (array1[i]->num == array2[i]->num && array1[i]->line == array2[i]->line) {
                trueCounter++;
            }
        }
        cout << "Equals object = " << trueCounter << endl;

        for (int i = 0; i < n; i++) {
            delete array1[i];
        }
        delete[] array1;
        delete[] array2;
        //delete bh;
    }

    return 0;
}