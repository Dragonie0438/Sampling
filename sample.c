#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <ctype.h>
#include <windows.h>
#include <psapi.h>

#define M_PI 3.14159265358979323846 
#define NUM_SAMPLES 1000000      // 采样总次数，可以自行调整
#define SQRT_E 1.6487212707     // e的平方根
#define KARNEY_EPSILON 1e-12           // 浮点误差阈值，消除因浮点运算精度不足而带来的误差


// ================== 内存监控（此部分在代码中作为内存检测作用，并非主要代码部分） ==================
size_t get_current_memory_kb() {
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        return pmc.WorkingSetSize / 1024;
    }
    return 0;
}

void print_memory_usage(const char* phase) {
    size_t mem_kb = get_current_memory_kb();
    printf(" [内存监控] %-20s: %6lu KB\n", phase, mem_kb);
}
// ================整数高斯分布函数===============

double normal_pdf(int k, double mu, double sigma ,double x) {
    double z = (k - mu) / sigma;
    double ρ1 = exp(-z * z / 2.0) / (sigma * sqrt(2.0 * M_PI));// 标准正态分布公式
    return ρ1 / x; //归一化
}

double sum_pdf(double mu, double sigma, int lower_bound, int upper_bound) {
    int i = 0;
    double sum = 0.0 , z;
    for (int i = lower_bound; i <= upper_bound; i++) {
        z = (i - mu) / sigma;
        sum += exp(-z * z / 2.0) / (sigma * sqrt(2.0 * M_PI));
    }
    return sum;
} //用于整数正态分布的归一化


// ================== 拒绝采样法 ==================

int sample_by_rejection(double mu, double sigma) {

    int lower_bound = (int)floor(mu - 5 * sigma);
    int upper_bound = (int)ceil(mu + 5 * sigma);  // 取整数范围
    int range = upper_bound - lower_bound + 1;
    double ρ2 = sum_pdf(mu, sigma, lower_bound, upper_bound) ;
    int m = (int)round(mu); // 概率最高点
    double max_p = normal_pdf(m, mu, sigma,ρ2);

   
    while (1) {
        int k = rand() % range + lower_bound;
        double u = (double)rand() / RAND_MAX;  // 随机生成一个[0,1]之间的随机数
        if (u < normal_pdf(k, mu, sigma,ρ2) / max_p) {
            return k;
        }  // 拒绝采样算法的核心部分，筛选出符合概率分布的样本
    }
}

// ================== 查表法 ==================
typedef struct {
    int* values; // 存储所有整数值
    double* cdf;  // 累积概率函数表的数组
    int size; // 查表容量
    int lower_bound; //最小整数值
} LookupTable;

LookupTable* init_lookup_table(double mu, double sigma, double threshold) {
    // 确定范围（±5σ）
    int lower_bound = (int)floor(mu - 5 * sigma);
    int upper_bound = (int)ceil(mu + 5 * sigma);
    double ρ2 = sum_pdf(mu, sigma, lower_bound, upper_bound) ;

    // 分配内存
    int size = upper_bound - lower_bound + 1;
    LookupTable* table = (LookupTable*)malloc(sizeof(LookupTable));
    table->values = (int*)malloc(size * sizeof(int));
    table->cdf = (double*)malloc(size * sizeof(double));
    table->size = size;
    table->lower_bound = lower_bound;

    // 构建累积分布表
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

    return table;
}

// 查表法采样
int sample_by_lookup(LookupTable* table) {
    double u = (double)rand() / RAND_MAX;

    // 二分查找定位样本
    int left = 0, right = table->size - 1;
    while (left < right) {
        int mid = left + (right - left) / 2;
        u < table->cdf[mid] ? (right = mid) : (left = mid + 1);
    }
    return table->values[left];
}

// 释放查表法内存
void free_lookup_table(LookupTable* table) {
    free(table->values);
    free(table->cdf);
    free(table);
}

// ================== Karney算法1（使用拒绝采样法进行k的获取） ==================
int Karney_sample_1(double mu, double sigma) {
    const int sigma_ceil = (int)ceil(sigma);
    const double sigma_inv = 1.0 / sigma;
    double p1, p2;
    int k;

    while (1) {

        while (1) {
            k = rand() % 5; 
            p1 = 0.393469340287366577 * exp(-0.5 * k);
            if ((double)rand() / (RAND_MAX + 1.0) < p1) {     //拒绝采样的方式获取k
                p2 = exp(-0.5 * k * (k - 1));
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
// ================== Karney算法2 （使用查表法进行k值的获取）==================
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
        double prob = C * exp(-lambda * k);// k获取的分布概率密度函数
        total += prob;
        et->cum_probs[i] = total;
    }

    for (int i = 0; i < et->table_size; ++i) {
        et->cum_probs[i] /= total;
    }

    return et;
}


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
//以上的查表方式与CDF查表法类似

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

// ================== 主函数 ==================
void clear_input_buffer() {
    int c;
    while ((c = getchar()) != '\n' && c != EOF);
}

int get_user_choice() {     //选择采样方法
    int choice = 0;
    printf(" 1.拒绝采样法\n");
    printf(" 2.查表法\n");
    printf(" 3.Karney算法-1\n");
    printf(" 4.Karney算法-2\n");
    printf(" 输入选择(1-4): ");

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
    srand((unsigned int)time(NULL));  // 初始化随机种子

    double mu, sigma;
    printf(" μ: ");
    scanf_s("%lf", &mu);
    printf(" σ: ");
    scanf_s("%lf", &sigma);
    clear_input_buffer();

    // 选择采样方法
    int method_choice = get_user_choice();
    LookupTable* table = NULL;  // 查表法专用指针

    // 根据选择初始化查表法
    if (method_choice == 2) {
        table = init_lookup_table(mu, sigma, 1e-6);
    }

    // 打开输出文件
    FILE* file = fopen("C:/Users/17839/Desktop/Sample/samples.txt", "w");// 输出文件路径，可自行修改
    if (!file) {
        printf("无法打开文件\n");
        if (table) free_lookup_table(table);
        return 1;
    }
    fprintf(file, "%.6f\n%.6f\n", mu, sigma);   // 将mu和sigma的参数也写入文件中，便于分析，可注释掉，不影响采样结果

    printf("\n 开始采样...\n");
    clock_t start = clock();


    for (int i = 0; i < NUM_SAMPLES; i++) {
        int sample;

        switch (method_choice) {
        case 1:
            sample = sample_by_rejection(mu, sigma);
            break;

        case 2:
            sample = sample_by_lookup(table);
            break;

        case 3:
            sample = Karney_sample_1(mu, sigma);
            break;

        case 4:
            sample = Karney_sample_2(mu, sigma);
            break;
        }

        fprintf(file, "%d\n", sample);

        if (i % 100000 == 0 && i > 0) {
            printf(" 已生成 %d 个样本...\n", i);
            if (i % 200000 == 0) print_memory_usage(" 采样中");
        }
    }

    clock_t end = clock();
    fclose(file);

    printf("\n 采样完成 \n");
    double duration = (double)(end - start) / CLOCKS_PER_SEC;
    printf(" 耗时: %.2f秒\n", duration);
    printf(" 速度: %.0f次/秒\n", NUM_SAMPLES / duration);

    if (table) free_lookup_table(table);
    return 0;
}