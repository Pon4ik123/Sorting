#include <iostream>
#include <math.h>

using namespace std;

template<typename T>
class SomeObject {
public:
    SomeObject() = default;
    SomeObject(T D, T C) : d(D), c(C) {}
    ~SomeObject() {}

    T d;
    T c;
};

template<typename T>
class MaxHeap {
public:
    T *array;
    int size;
    int capacity;

    MaxHeap() : size(0), capacity(1) {
        array = new T[capacity];
    }

    MaxHeap(T array[], int capacity) {
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

    void heapify(int index) {
        int l = left(index);
        int r = right(index);
        int max = index;

        if (r < size && array[r] > array[max]) {
            max = r;
        }
        if (l < size && array[l] > array[max]) {
            max = l;
        }
        if (max != index) {
            swap(array[index], array[max]);
            heapify(max);
        }
    }

    void sort(T* arr, int n) {
        for (int i = n / 2 - 1; i >= 0; i--){
            heapify(arr, n, i);
        }

        for (int i = n - 1; i >= 0; i--) {
            swap(arr[0], arr[i]);
            heapify(arr, i, 0);
        }
    }

    void deleteHeap() {
        delete[] array;
        array = nullptr;
        size = 0;
        capacity = 0;
    }
};

void countingSortInteger(int *array, int size, int m) {
    int *count = new int[m]();
    int *result = new int[size];

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

template<typename T>
void bubbleSort(T *array, int size) {
    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - i - 1; j++) {
            if (array[j] > array[j + 1]) {
                T temp = array[j];
                array[j] = array[j + 1];
                array[j + 1] = temp;
            }
        }
    }
}

void bucketSortInt(int *array, int n, int m) {
    if (array == nullptr || n <= 1) {
        return;
    }

    int **buckets = new int *[m];
    for (int i = 0; i < m; ++i) {
        buckets[i] = new int[n];
    }

    int *bucketSizes = new int[m]();

    for (int i = 0; i < n; ++i) {
        int index = array[i] * m / (std::numeric_limits<int>::max() + 1);
        buckets[index][bucketSizes[index]++] = array[i];
    }

    for (int i = 0; i < m; ++i) {
        bubbleSort(buckets[i], *(buckets[i] + bucketSizes[i]));
    }

    int index = 0;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < bucketSizes[i]; ++j) {
            array[index++] = buckets[i][j];
        }
    }

    for (int i = 0; m; ++i) {
        delete[] buckets[i];
    }
    delete[] buckets;
    delete[] bucketSizes;
}

template<typename T>
void bucketSortObject(SomeObject<T> *array, int n, int m, T(*keyFunc)(SomeObject<T>*), bool(*comparator)(SomeObject<T>*, SomeObject<T>*)) {
    if (array == nullptr || n <= 1 || keyFunc == nullptr || comparator == nullptr) {
        return;
    }

    SomeObject<T> **buckets = new SomeObject<T> *[m];
    for (int i = 0; i < m; ++i) {
        buckets[i] = new SomeObject<T>[n];
    }

    int *bucketSizes = new int[m]();

    for (int i = 0; i < n; ++i) {
        T key = keyFunc(&array[i]);
        int index = key * m / (std::numeric_limits<T>::max() + 1);
        buckets[index][bucketSizes[index]++] = array[i];
    }

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < bucketSizes[i]; ++j) {
            for (int k = j + 1; k < bucketSizes[i]; ++k) {
                if (comparator(&buckets[i][k], &buckets[i][j])) {
                    bubbleSort(buckets[i][k], buckets[i][j]);
                }
            }
        }
    }

    int index = 0;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < bucketSizes[i]; ++j) {
            array[index++] = buckets[i][j];
        }
    }

    for (int i = 0; i < m; ++i) {
        delete[] buckets[i];
    }
    delete[] buckets;
    delete[] bucketSizes;
}

template<typename T>
void printArray(T array[], int size) {
    for (int i = 0; i < size; i++) {
        cout << array[i] << " ";
    }
    cout << endl;
}

int main() {
    return 0;
}
