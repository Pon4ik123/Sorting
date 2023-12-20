#include <iostream>
#include <algorithm>

using namespace std;

template <typename T>
class MaxHeap {
public:
    T* array;
    int size;
    int capacity;

    MaxHeap() : size(0), capacity(1) {
        array = new T[capacity];
    }

    MaxHeap(T array[], int capacity){
        this->array = new T[capacity];
        this->capacity = capacity;
        for (int i = 0; i < capacity; ++i) {
            add(array[i]);
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

    void swap(T &x, T &y) {
        T temp = x;
        x = y;
        y = temp;
    }

    void add(T value) {
        if (size == capacity) {
            cout << "Heap is full" << endl;
            return;
        }
        array[size] = value;
        heapifyUp(size);
        size++;
    }


    void heapifyUp(int index) {
        while (index != 0 && array[parent(index)] < array[index]) {
            swap(array[index], array[parent(index)]);
            index = parent(index);
        }
    }

    void heapifyDown(int index) {
        int l = left(index);
        int r = right(index);
        int max = index;

        if (r<size && array[r]>array[max]) {
            max = r;
        }
        if (l < size && array[l] > array[max]) {
            max = l;
        }
        if (max != index) {
            swap(array[index], array[max]);
            heapifyDown(max);
        }
    }

    void deleteHeap() {
        delete[] array;
        array = nullptr;
        size = 0;
        capacity = 0;
    }
};

template <typename T>
void printArray(T array[], int size) {
    for (int i = 0; i < size; i++)
        cout << array[i] << " ";
    cout << endl;
}

template <typename T>
void countingSort(T array[], int size){
    int b[size];
    int max = array[0];
    for(int i = 1; i < size; i++){
        if(array[i]>max){
            max = array[i];
        }
    }

    int count[max+1];
    for(int i = 0; i <= max; i++){
        count[i] = 0;
    }

    for(int i = 0; i < size; i++){
        count[array[i]]++;
    }

    for(int i = 1; i <= max; i++){
        count[i] += count[i-1];
    }

    for(int i = size-1; i >= 0; i--){
        b[--count[array[i]]] = array[i];
    }

    for(int i = 0; i < size; i++){
        array[i] = b[i];
    }
}

template <typename T>
void bubleSort(T array[], int size){
    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - i - 1; j++) {
            if (array[j] > array[j + 1]) {
                T temp = array[j];
                array[j] = array[j + 1];
                array[j + 1] = temp;
            }
        }
    }
    cout << "Table has been sorted using bubble sort." << endl;
}

template <typename T>
void bucketSort(T array[], int size){
    int maxVal = array[0];
    for (int i = 1; i < size; i++) {
        if (array[i] > maxVal) {
            maxVal = array[i];
        }
    }

    // Create n empty buckets
    const int numBuckets = size;
    int* bucketSizes = new int[numBuckets](); // Array to store the size of each bucket
    int** buckets = new int*[numBuckets]; // Array of arrays to represent buckets

    for (int i = 0; i < numBuckets; i++) {
        buckets[i] = new int[size];
    }

    // Put array elements in different buckets
    for (int i = 0; i < size; i++) {
        int bi = static_cast<int>(size * (static_cast<float>(array[i]) / maxVal));
        buckets[bi][bucketSizes[bi]++] = array[i];
    }

    // Sort individual buckets using std::sort
    for (int i = 0; i < numBuckets; i++) {
        bubleSort(buckets[i], bucketSizes[i]);
    }

    // Concatenate all buckets into arr[]
    int index = 0;
    for (int i = 0; i < numBuckets; i++) {
        for (int j = 0; j < bucketSizes[i]; j++) {
            array[index++] = buckets[i][j];
        }
    }

    // Clean up dynamically allocated memory
    delete[] bucketSizes;
    for (int i = 0; i < numBuckets; i++) {
        delete[] buckets[i];
    }
    delete[] buckets;
}