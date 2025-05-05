#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <ctype.h>
#include <windows.h>
#include <psapi.h>

#define M_PI 3.14159265358979323846
#define NUM_SAMPLES 1000000
#define SQRT_E 1.6487212707     // sqrt(e)
#define KARNEY_EPSILON 1e-12           // 浮点误差阈值

// ================== 内存监控函数 ==================
size_t get_current_memory_kb() {
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        return pmc.WorkingSetSize / 1024;
    }
    return 0;
}

void print_memory_usage(const char* phase) {
    size_t mem_kb = get_current_memory_kb();
    printf("[内存监控] %-20s: %6lu KB\n", phase, mem_kb);
}

// ================== 拒绝采样法 ==================
double normal_pdf(int k, double mu, double sigma ,double x) {
    double z = (k - mu) / sigma;
    double ρ1 = exp(-z * z / 2.0) / (sigma * sqrt(2.0 * M_PI));
    return ρ1 / x;
}


double sum_pdf(double mu, double sigma, int lower_bound, int upper_bound) {
    int i = 0;
    double sum = 0.0 , z;
    for (int i = lower_bound; i <= upper_bound; i++) {
        z = (i - mu) / sigma;
        sum += exp(-z * z / 2.0) / (sigma * sqrt(2.0 * M_PI));
    }
    return sum;
}

int sample_by_rejection(double mu, double sigma) {

    int lower_bound = (int)floor(mu - 5 * sigma);
    int upper_bound = (int)ceil(mu + 5 * sigma);
    int range = upper_bound - lower_bound + 1;
    double ρ2 = sum_pdf(mu, sigma, lower_bound, upper_bound) ;
    int mode = (int)round(mu);
    double max_p = normal_pdf(mode, mu, sigma,ρ2);

    while (1) {
        int k = rand() % range + lower_bound;
        double u = (double)rand() / RAND_MAX;
        if (u < normal_pdf(k, mu, sigma,ρ2) / max_p) {
            return k;
        }
    }
}

// ================== 查表法 ==================
typedef struct {
    int* values;
    double* cdf;
    int size;
    int lower_bound;
} LookupTable;

LookupTable* init_lookup_table(double mu, double sigma, double threshold) {
    print_memory_usage("查表初始化前");
    int lower_bound = (int)floor(mu - 5 * sigma);
    int upper_bound = (int)ceil(mu + 5 * sigma);

    double ρ2 = sum_pdf(mu, sigma, lower_bound, upper_bound) ;


    while (normal_pdf(lower_bound, mu, sigma, ρ2) < threshold) lower_bound++;
    while (normal_pdf(upper_bound, mu, sigma , ρ2) < threshold) upper_bound--;

    int size = upper_bound - lower_bound + 1;
    printf("查表范围: %d 到 %d (大小: %d)\n", lower_bound, upper_bound, size);

    LookupTable* table = (LookupTable*)malloc(sizeof(LookupTable));
    table->values = (int*)malloc(size * sizeof(int));
    table->cdf = (double*)malloc(size * sizeof(double));
    table->size = size;
    table->lower_bound = lower_bound;

    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        int k = lower_bound + i;
        table->values[i] = k;
        sum += normal_pdf(k, mu, sigma, ρ2);
        table->cdf[i] = sum;
    }

    for (int i = 0; i < size; i++) {
        table->cdf[i] /= sum;
    }

    print_memory_usage("查表初始化后");
    return table;
}

int sample_by_lookup(LookupTable* table) {

    double u = (double)rand() / RAND_MAX;
    int left = 0, right = table->size - 1;
    while (left < right) {
        int mid = left + (right - left) / 2;
        u < table->cdf[mid] ? (right = mid) : (left = mid + 1);
    }
    return table->values[left];
}

void free_lookup_table(LookupTable* table) {
    free(table->values);
    free(table->cdf);
    free(table);
}

// ================== Karney算法1 ==================

int Karney_sample_1(double mu, double sigma) {
    const int sigma_ceil = (int)ceil(sigma);
    const double sigma_inv = 1.0 / sigma;
    double p1, p2;
    int k;

    while (1) {

        while (1) {
            k = rand() % 5; // 选择整数k >= 0
            p1 = 0.393469340287366577 * exp(-0.5 * k);// 计算选择k的概率p1
            if ((double)rand() / (RAND_MAX + 1.0) < p1) {
                p2 = exp(-0.5 * k * (k - 1));// 计算接受k的概率p2
                if ((double)rand() / (RAND_MAX + 1.0) < p2) {
                    break;
                }
            }
        }
    
        int s = (rand() % 2) ? 1 : -1;

        double shifted = sigma * k + s * mu;
        int i0 = (int)floor(shifted + KARNEY_EPSILON);
        double x0 = (i0 - shifted) * sigma_inv;
        int j = rand() % sigma_ceil;
        double x = x0 + j * sigma_inv;

        if (x >= 1.0 - KARNEY_EPSILON) continue;
        if (k == 0 && fabs(x) < KARNEY_EPSILON && s < 0) continue;
        if ((double)rand() / RAND_MAX > exp(-0.5 * x * (2 * k + x))) continue;

        return s * (i0 + j);
    }
}
// ================== Karney算法2 ==================
typedef struct {
    int* k_values;      // 存储离散值k
    double* cum_probs;  // 累积概率数组
    int table_size;     // 查表容量
    int base_k;         // 最小k值
} ExpLookupTable;

// 初始化指数分布查表（k=0到6）
ExpLookupTable* init_exp_lookup() {
    const double C = 0.393469340287366577;
    const double lambda = 0.5;
    const int min_k = 0;
    const int max_k = 6;

    ExpLookupTable* et = malloc(sizeof(ExpLookupTable));
    et->base_k = min_k;
    et->table_size = max_k - min_k + 1;

    et->k_values = malloc(et->table_size * sizeof(int));
    et->cum_probs = malloc(et->table_size * sizeof(double));

    double total = 0.0;
    for (int i = 0; i < et->table_size; ++i) {
        int k = min_k + i;
        et->k_values[i] = k;
        double prob = C * exp(-lambda * k);
        total += prob;
        et->cum_probs[i] = total;
    }

    for (int i = 0; i < et->table_size; ++i) {
        et->cum_probs[i] /= total;
    }

    return et;
}

// 指数分布离散采样
int exp_discrete_sample(ExpLookupTable* et) {
    double u = (double)rand() / RAND_MAX;
    int low = 0, high = et->table_size - 1;
    while (low < high) {
        int mid = low + (high - low) / 2;
        if (u < et->cum_probs[mid]) {
            high = mid;
        }
        else {
            low = mid + 1;
        }
    }
    return et->k_values[low];
}

void free_exp_lookup(ExpLookupTable* et) {
    if (et) {
        free(et->k_values);
        free(et->cum_probs);
        free(et);
    }
}

int Karney_sample_2(double mu, double sigma) {
    static ExpLookupTable* et = NULL;
    if (!et) et = init_exp_lookup();  // 单次初始化

    const int sigma_ceil = (int)ceil(sigma);
    const double sigma_inv = 1.0 / sigma;

    while (1) {
        
        int k = exp_discrete_sample(et);
        double p2 = exp(-0.5 * k * (k - 1));

        if ((double)rand() / (RAND_MAX + 1.0) >= p2) continue;

        int s = (rand() % 2) ? 1 : -1;

        double shifted = sigma * k + s * mu;
        int i0 = (int)floor(shifted + KARNEY_EPSILON);
        double x0 = (i0 - shifted) * sigma_inv;

        int j = rand() % sigma_ceil;
        double x = x0 + j * sigma_inv;

        if (x >= 1.0 - KARNEY_EPSILON) continue;
        if (k == 0 && fabs(x) < KARNEY_EPSILON && s < 0) continue;
        if ((double)rand() / RAND_MAX > exp(-0.5 * x * (2 * k + x))) continue;

        return s * (i0 + j);
    }
}

// ================== 用户界面 ==================
void clear_input_buffer() {
    int c;
    while ((c = getchar()) != '\n' && c != EOF);
}

int get_user_choice() {
    int choice = 0;
    printf("\n请选择采样算法:\n");
    printf("1. 拒绝采样法\n");
    printf("2. 查表法\n");
    printf("3. Karney算法1\n");
    printf("4. Karney算法2\n");
    printf("请输入选择(1-4): ");

    while (1) {
        if (scanf_s("%d", &choice) && choice >= 1 && choice <= 4) {
            clear_input_buffer();
            return choice;
        }
        clear_input_buffer();
        printf("无效输入，请重新输入(1-4): ");
    }
}

int main() {
    srand((unsigned int)time(NULL));
    print_memory_usage("程序启动后");

    const char* output_path = "C:/Users/17839/Desktop/Sample/samples.txt";
    printf("整数高斯采样程序 (固定采样100万次)\n");
    printf("--------------------------------\n");

    double mu, sigma;
    printf(".请输入均值μ: ");
    scanf_s("%lf", &mu);
    printf("请输入标准差σ: ");
    scanf_s("%lf", &sigma);
    clear_input_buffer();

    int method_choice = get_user_choice();
    LookupTable* table = NULL;
    bool use_lookup = (method_choice == 2);
    bool use_karney_1 = (method_choice == 3);
    bool use_karney_2 = (method_choice == 4);

    if (use_lookup) {
        table = init_lookup_table(mu, sigma, 1e-6);
    }

    FILE* file = fopen(output_path, "w");
    if (!file) {
        printf("无法打开文件\n");
        if (table) free_lookup_table(table);
        return 1;
    }
    fprintf(file, "%.6f\n%.6f\n", mu, sigma);

    printf("\n开始采样...\n");
    clock_t start = clock();
    print_memory_usage("采样开始前");

    for (int i = 0; i < NUM_SAMPLES; i++) {
        int sample = use_lookup ? sample_by_lookup(table) :
            use_karney_1 ? Karney_sample_1(mu, sigma) :
            use_karney_2 ? Karney_sample_2(mu, sigma) :
            sample_by_rejection(mu, sigma);
        
        fprintf(file, "%d\n", sample);

        if (i % 100000 == 0 && i > 0) {
            printf("已生成 %d 个样本...\n", i);
            if (i % 200000 == 0) print_memory_usage("采样中");
        }
    }

    clock_t end = clock();
    fclose(file);

    printf("\n采样完成!\n");
    printf("耗时: %.2f秒  速度: %.0f次/秒\n",
        (double)(end - start) / CLOCKS_PER_SEC,
        NUM_SAMPLES / ((double)(end - start) / CLOCKS_PER_SEC));

    if (table) {
        free_lookup_table(table);
        print_memory_usage("释放查表后");
    }
    print_memory_usage("程序结束前");
    return 0;
}