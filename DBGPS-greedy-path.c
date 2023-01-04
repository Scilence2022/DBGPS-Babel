#include <stdio.h>
//#include <stdlib.h>
//#include <conio.h>
//#include <curses.h>
#include <stdint.h>
#include <zlib.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>
#include "ketopt.h" // command-line argument parser
#include "kthread.h" // multi-threading models: pipeline and multi-threaded for loop


#include "kseq.h" // FASTA/Q parser
KSEQ_INIT(gzFile, gzread)

#include "khashl.h" // hash table
#define KC_BITS 10
#define KC_MAX ((1<<KC_BITS) - 1)
#define kc_c4_eq(a, b) ((a)>>KC_BITS == (b)>>KC_BITS) // lower 10 bits for counts; higher bits for k-mer
#define kc_c4_hash(a) ((a)>>KC_BITS)

#define Max_Path_Num 100
#define Max_Path_Len 300
#define Max_Cov_Ratio  5


KHASHL_SET_INIT(, kc_c4_t, kc_c4, uint64_t, kc_c4_hash, kc_c4_eq)

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

const unsigned char nt4_seq_table[4] = { // translate 0123 to ACGT
        'A', 'C', 'G', 'T'};


static inline uint64_t hash64(uint64_t key, uint64_t mask) // invertible integer hash function
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

// The inversion of hash64(). Modified from <https://naml.us/blog/tag/invertible>
static inline uint64_t hash64i(uint64_t key, uint64_t mask) //
{
    uint64_t tmp;
    // Invert key = key + (key << 31)
    tmp = (key - (key << 31)); 	key = (key - (tmp << 31)) & mask;
    // Invert key = key ^ (key >> 28)
    tmp = key ^ key >> 28; 	key = key ^ tmp >> 28;
    // Invert key *= 21
    key = (key * 14933078535860113213ull) & mask;
    // Invert key = key ^ (key >> 14)
    tmp = key ^ key >> 14; 	tmp = key ^ tmp >> 14; 	tmp = key ^ tmp >> 14; 	key = key ^ tmp >> 14;
    // Invert key *= 265
    key = (key * 15244667743933553977ull) & mask;
    // Invert key = key ^ (key >> 24)
    tmp = key ^ key >> 24; 	key = key ^ tmp >> 24;
    // Invert key = (~key) + (key << 21)
    tmp = ~key; 	tmp = ~(key - (tmp << 21)); tmp = ~(key - (tmp << 21)); key = ~(key - (tmp << 21)) & mask;
    return key;
}


unsigned char* uint64_acgt(uint64_t key, unsigned char* seq, unsigned char km_len){ //decode uint_64_t kmer  to actg
    int p = km_len-1;
    while (p >= 0){
        int n = key % 4;
        key = key >> 2;
        *(seq + p) = nt4_seq_table[n];
        p = p - 1;
    }
    return seq;
}

unsigned char* uint64_int8(uint64_t key, unsigned char* seq){
    //decode uint_64_t kmer  to  char

    uint64_t mask = (1ULL<<8)-1;
    int i =0;
    for(i =0; i< 8; i++){
        *(seq+7-i) = key & mask;
        mask = mask <<8;
    }
    return seq;
}


const unsigned char two_to_one_table[256] = { //one base encoder to two base encoder
    51, 50, 49, 48, 60, 61, 62, 63, 58, 56, 59, 57, 53, 55, 52, 54,
    35, 34, 33, 32, 44, 45, 46, 47, 42, 40, 43, 41, 37, 39, 36, 38,
    19, 18, 17, 16, 28, 29, 30, 31, 26, 24, 27, 25, 21, 23, 20, 22,
    3, 2, 1, 0, 12, 13, 14, 15, 10, 8, 11, 9, 5, 7, 4, 6,
    195, 194, 193, 192, 204, 205, 206, 207, 202, 200, 203, 201, 197, 199, 196, 198,
    211, 210, 209, 208, 220, 221, 222, 223, 218, 216, 219, 217, 213, 215, 212, 214,
    227, 226, 225, 224, 236, 237, 238, 239, 234, 232, 235, 233, 229, 231, 228, 230,
    243, 242, 241, 240, 252, 253, 254, 255, 250, 248, 251, 249, 245, 247, 244, 246,
    163, 162, 161, 160, 172, 173, 174, 175, 170, 168, 171, 169, 165, 167, 164, 166,
    131, 130, 129, 128, 140, 141, 142, 143, 138, 136, 139, 137, 133, 135, 132, 134,
    179, 178, 177, 176, 188, 189, 190, 191, 186, 184, 187, 185, 181, 183, 180, 182,
    147, 146, 145, 144, 156, 157, 158, 159, 154, 152, 155, 153, 149, 151, 148, 150,
    83, 82, 81, 80, 92, 93, 94, 95, 90, 88, 91, 89, 85, 87, 84, 86,
    115, 114, 113, 112, 124, 125, 126, 127, 122, 120, 123, 121, 117, 119, 116, 118,
    67, 66, 65, 64, 76, 77, 78, 79, 74, 72, 75, 73, 69, 71, 68, 70,
    99, 98, 97, 96, 108, 109, 110, 111, 106, 104, 107, 105, 101, 103, 100, 102
};


const unsigned char one_to_two_table[256] = { //one base encoder to two base encoder

    51, 50, 49, 48, 62, 60, 63, 61, 57, 59, 56, 58, 52, 53, 54, 55,
    35, 34, 33, 32, 46, 44, 47, 45, 41, 43, 40, 42, 36, 37, 38, 39,
    19, 18, 17, 16, 30, 28, 31, 29, 25, 27, 24, 26, 20, 21, 22, 23,
    3, 2, 1, 0, 14, 12, 15, 13, 9, 11, 8, 10, 4, 5, 6, 7,
    227, 226, 225, 224, 238, 236, 239, 237, 233, 235, 232, 234, 228, 229, 230, 231,
    195, 194, 193, 192, 206, 204, 207, 205, 201, 203, 200, 202, 196, 197, 198, 199,
    243, 242, 241, 240, 254, 252, 255, 253, 249, 251, 248, 250, 244, 245, 246, 247,
    211, 210, 209, 208, 222, 220, 223, 221, 217, 219, 216, 218, 212, 213, 214, 215,
    147, 146, 145, 144, 158, 156, 159, 157, 153, 155, 152, 154, 148, 149, 150, 151,
    179, 178, 177, 176, 190, 188, 191, 189, 185, 187, 184, 186, 180, 181, 182, 183,
    131, 130, 129, 128, 142, 140, 143, 141, 137, 139, 136, 138, 132, 133, 134, 135,
    163, 162, 161, 160, 174, 172, 175, 173, 169, 171, 168, 170, 164, 165, 166, 167,
    67, 66, 65, 64, 78, 76, 79, 77, 73, 75, 72, 74, 68, 69, 70, 71,
    83, 82, 81, 80, 94, 92, 95, 93, 89, 91, 88, 90, 84, 85, 86, 87,
    99, 98, 97, 96, 110, 108, 111, 109, 105, 107, 104, 106, 100, 101, 102, 103,
    115, 114, 113, 112, 126, 124, 127, 125, 121, 123, 120, 122, 116, 117, 118, 119


};


 unsigned char* one_encode_two(unsigned char* seq, int seq_len){

    int i;
    for(i=0; i < seq_len; i++){

        seq[i] = one_to_two_table[i];
    }
    return seq;
}


unsigned char* two_encode_one(unsigned char* seq, int seq_len){

    int i;
    for(i=0; i < seq_len; i++){
        seq[i] = two_to_one_table[i];
    }
    return seq;
}


uint64_t comp_rev2(uint64_t x, unsigned char km_len){
    uint64_t a, b=0ULL, mask = (1ULL<<km_len*2) - 1;

    a = ~x; a = a & mask;
    while((a & mask) > 0 ){
        b = b<<2;
        b = b + (a & 3ULL);
        a = a>>2;
    }
//    b = ~b;
    b = b & mask;
    return b;
}


uint64_t comp_rev(uint64_t x, unsigned char km_len){
//    if (c < 4) { // not an "N" base
//        x[0] = (x[0] << 2 | c) & mask;                  // forward strand
//        x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
//    uint64_t y=0;// mask = (1ULL << km_len * 2) - 1;
    uint64_t y=0;
    int i, c;
    for(i=0;i<km_len;i++){
        y = y << 2;
        c = x & 3ULL;
        y = y | (3-c);
        x = x >> 2;
    }
    return y;
}


uint64_t min_hash_key(uint64_t x, unsigned char km_len){
    uint64_t y = comp_rev(x, km_len);
    return (x < y) ? x:y;
}


uint64_t actgkmer_uint64(unsigned char* seq, unsigned char km_len){
//    Just for testing
//// Function verified 20210708 Lifu Song
    int i;
    uint64_t x, mask = (1ULL<<km_len*2) - 1; //shift = (km_len - 1) * 2
    for (i = 0, x = 0; i < km_len; ++i) {
        //fprintf(stdout, "%c ", seq[i]);
        int c = seq_nt4_table[(uint8_t)seq[i]];
        //fprintf(stdout, "%d ", c);
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & mask;                  // forward strand
//            x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
        } else x = 0; // if there is an "N", restart
    }

    //fprintf(stdout, "y %ld\n", x);
    return x;
}


uint64_t actgkmer_hashkey(unsigned char* seq, unsigned char km_len){
//// Function verified 20210708 Lifu Song
    int i;
    uint64_t x[2], y,  shift = (km_len - 1) * 2, mask = (1ULL<<km_len*2) - 1;
    for (i = 0, x[0] = x[1] = 0; i < km_len; ++i) {
        //fprintf(stdout, "%c ", seq[i]);
        int c = seq_nt4_table[(uint8_t)seq[i]];
        //fprintf(stdout, "%d ", c);
        if (c < 4) { // not an "N" base
            x[0] = (x[0] << 2 | c) & mask;                  // forward strand
            x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
        } else x[0] = x[1] = 0; // if there is an "N", restart
    }
    y = x[0] < x[1]? x[0] : x[1];
    //fprintf(stdout, "y %ld\n", y);
    return y;
}

unsigned long uint8_actg(unsigned char bts){
    unsigned int i=0;
    unsigned long four_base_seq=0;

    for(i =0 ; i < 4; ++i){
        int m = bts & 3;
        bts = bts >>2;
        four_base_seq =four_base_seq >>8;
        four_base_seq = four_base_seq + ((unsigned long)nt4_seq_table[m]<< 24);
    }
    return four_base_seq;
}


unsigned char actg_uint8(unsigned long actg_seq){
    unsigned int i=0;
    unsigned char uint8_value=0;

    for(i =0 ; i < 4; ++i){
        unsigned int m = actg_seq & 255;
        actg_seq = actg_seq >> 8;
        uint8_value = uint8_value >> 2;
        uint8_value = uint8_value + ((unsigned char)seq_nt4_table[m]<<6);
    }
    return uint8_value;
}

void actg_intseq(unsigned char* actg_seq, unsigned char* int_seq, int8_t actg_seq_len){

    int p = 0;
    unsigned long n, n1, n2, n3, n4;
    while(p < actg_seq_len-3){
        n1 = actg_seq[p];
        n2 = actg_seq[p+1];
        n3 = actg_seq[p+2];
        n4 = actg_seq[p+3];
        n = (n1<<24) + (n2<<16) + (n3<<8) + n4;
        int_seq[p>>2] = actg_uint8(n) ;
        p = p + 4;
    }
    int_seq[p + 1] = '\0';
}


void intseq_actg(unsigned char* int_seq, unsigned char* actg_seq, int8_t int_seq_len){
//// int seq to DNA seq
//// Function not tested yet
//// one char(byte) to four bases
    int p = 0;
    unsigned long n, mask;

    unsigned char n1, n2, n3, n4;
    while(p < int_seq_len){
        mask= (1<<8) - 1;
        n = uint8_actg(int_seq[p]) ;
        n4 = (n & mask);
        mask = mask <<8;
        n3 = (n & mask)>>8;
        mask = mask <<8;
        n2 = (n & mask)>>16;
        mask = mask <<8;
        n1 = (n & mask)>>24;

        actg_seq[p<<2] = n1;
        actg_seq[(p<<2) + 1] = n2;
        actg_seq[(p<<2) + 2] = n3;
        actg_seq[(p<<2) + 3] = n4;
        p = p + 1;
    }
}

typedef struct {
    int p; // suffix length; at least 8
    kc_c4_t **h; // 1<<p hash tables
} kc_c4x_t;

static kc_c4x_t *c4x_init(int p)
{
    int i;
    kc_c4x_t *h;
    CALLOC(h, 1);
    MALLOC(h->h, 1<<p);
    h->p = p;
    for (i = 0; i < 1<<p; ++i)
        h->h[i] = kc_c4_init();
    return h;
}

typedef struct {
    int n, m;
    uint64_t *a;
} buf_c4_t;


static inline void c4x_insert_buf(buf_c4_t *buf, int p, uint64_t y) // insert a k-mer $y to a linear buffer
{
    int pre = y & ((1<<p) - 1);
    buf_c4_t *b = &buf[pre];
    if (b->n == b->m) {
        b->m = b->m < 8? 8 : b->m + (b->m>>1);
        REALLOC(b->a, b->m);
    }
    b->a[b->n++] = y;
}

static void count_seq_buf(buf_c4_t *buf, int k, int p, int len, const char *seq) // insert k-mers in $seq to linear buffer $buf
{
    int i, l;
    uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
    for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
        int c = seq_nt4_table[(uint8_t)seq[i]];
        if (c < 4) { // not an "N" base
            x[0] = (x[0] << 2 | c) & mask;                  // forward strand
            x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
            if (++l >= k) { // we find a k-mer
                uint64_t y = x[0] < x[1]? x[0] : x[1];
                c4x_insert_buf(buf, p, hash64(y, mask));
            }
        } else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
    }
}

typedef struct { // global data structure for kt_pipeline()
    int k, block_len, n_thread;
    kseq_t *ks;
    kc_c4x_t *h;
} pldat_t;

typedef struct { // data structure for each step in kt_pipeline()
    pldat_t *p;
    int n, m, sum_len, nk;
    int *len;
    char **seq;
    buf_c4_t *buf;
} stepdat_t;

static void worker_for(void *data, long i, int tid) // callback for kt_for()
{
    stepdat_t *s = (stepdat_t*)data;
    buf_c4_t *b = &s->buf[i];
    kc_c4_t *h = s->p->h->h[i];

    int j, p = s->p->h->p;
    for (j = 0; j < b->n; ++j) {
        khint_t k;
        int absent;
        k = kc_c4_put(h, b->a[j]>>p<<KC_BITS, &absent);
        if ((kh_key(h, k)&KC_MAX) < KC_MAX) ++kh_key(h, k);
    }
}

static void *worker_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    pldat_t *p = (pldat_t*)data;
    if (step == 0) { // step 1: read a block of sequences
        int ret;
        stepdat_t *s;
        CALLOC(s, 1);
        s->p = p;
        while ((ret = kseq_read(p->ks)) >= 0) {
            int l = p->ks->seq.l;
            if (l < p->k) continue;
            if (s->n == s->m) {
                s->m = s->m < 16? 16 : s->m + (s->n>>1);
                REALLOC(s->len, s->m);
                REALLOC(s->seq, s->m);
            }
            MALLOC(s->seq[s->n], l);
            memcpy(s->seq[s->n], p->ks->seq.s, l);
            s->len[s->n++] = l;
            s->sum_len += l;
            s->nk += l - p->k + 1;
            if (s->sum_len >= p->block_len)
                break;
        }
        if (s->sum_len == 0) free(s);
        else return s;
    } else if (step == 1) { // step 2: extract k-mers
        stepdat_t *s = (stepdat_t*)in;
        int i, n = 1<<p->h->p, m;
        CALLOC(s->buf, n);
        m = (int)(s->nk * 1.2 / n) + 1;
        for (i = 0; i < n; ++i) {
            s->buf[i].m = m;
            MALLOC(s->buf[i].a, m);
        }
        for (i = 0; i < s->n; ++i) {
            count_seq_buf(s->buf, p->k, p->h->p, s->len[i], s->seq[i]);
            free(s->seq[i]);
        }
        free(s->seq); free(s->len);
        return s;
    } else if (step == 2) { // step 3: insert k-mers to hash table
        stepdat_t *s = (stepdat_t*)in;
        int i, n = 1<<p->h->p;
        kt_for(p->n_thread, worker_for, s, n);
        for (i = 0; i < n; ++i) free(s->buf[i].a);
        free(s->buf); free(s);
    }
    return 0;
}

static kc_c4x_t *count_file(const char *fn, int k, int p, int block_size, int n_thread)
{
    pldat_t pl;
    gzFile fp;
    if ((fp = gzopen(fn, "r")) == 0) return 0;
    pl.ks = kseq_init(fp);
    pl.k = k;
    pl.n_thread = n_thread;
    pl.h = c4x_init(p);
    pl.block_len = block_size;
    kt_pipeline(3, worker_pipeline, &pl, 3);
    kseq_destroy(pl.ks);
    gzclose(fp);
    return pl.h;
}

static kc_c4x_t *count_file2(const char *fn, void *hh, int k, int p, int block_size, int n_thread)
{
    pldat_t pl;
    gzFile fp;
    if ((fp = gzopen(fn, "r")) == 0) return 0;
    pl.ks = kseq_init(fp);
    pl.k = k;
    pl.n_thread = n_thread;
    pl.h = (kc_c4x_t *)hh;
    pl.block_len = block_size;
    kt_pipeline(3, worker_pipeline, &pl, 3);
    kseq_destroy(pl.ks);
    gzclose(fp);
    return pl.h;
}


typedef struct {
    uint64_t c[256];
} buf_cnt_t;

typedef struct {
    const kc_c4x_t *h;
    buf_cnt_t *cnt;
} hist_aux_t;





unsigned short kmer_cov(uint64_t kmer, uint64_t mask, kc_c4x_t *h){

    int j, x, cov=0;
    uint64_t hash_key = hash64(kmer, mask);
    j = hash_key & ((1<<KC_BITS) - 1);
    if(kh_size(h->h[j]) < 1){return 0;}

    hash_key = hash_key >> KC_BITS<< KC_BITS;
    x = kc_c4_get(h->h[j], hash_key);
    //if(x == 0){return 50;} //To avoid of segment fault
    if(kh_exist(h->h[j], x)){
        cov = kh_key(h->h[j], x) & KC_MAX;
    }
    return cov;
}

typedef struct {
    struct path_node* pre_node;
    uint64_t kmer;
    unsigned int cov;
}path_node;

const uint64_t BI_MASK = 3ULL;


unsigned char four_kmers_to_byte(uint64_t k1, uint64_t k2, uint64_t k3, uint64_t k4){
    unsigned char m;
    m = (unsigned char)(((k1 & BI_MASK)<<6) + ((k2 & BI_MASK)<<4) + ((k3 & BI_MASK)<<2) + ((k4 & BI_MASK)));
    return m;
}



const uint64_t BYT_MASK = 255ULL;
unsigned char* nodes_path(path_node *term_node, unsigned char *path_seq, int path_len, int km_len, int index_byte_num){

    path_node *p_node;
    p_node = term_node;
    int path_seq_p = (path_len>>2) -1 + index_byte_num, nd_p=0;
    uint64_t kms[4];
    kms[nd_p] = p_node->kmer;

    int node_num=0;

    while(p_node->pre_node && path_seq_p >= index_byte_num){
        node_num++;
        p_node = p_node->pre_node;
        kms[nd_p] = p_node->kmer;

        if(nd_p >= 3){
            *(path_seq + path_seq_p) = four_kmers_to_byte(kms[3], kms[2], kms[1], kms[0]);
            nd_p = 0;
            path_seq_p--;
        }else{nd_p++;}
    }

    //putting index into path_seq

    while(p_node->pre_node){
        p_node = p_node->pre_node;
    }
    uint64_t index_kmer = p_node->kmer;
    int i;
    for(i=index_byte_num-1; i>=0; i--){
        *(path_seq + i) = (unsigned char) (index_kmer & BYT_MASK);
        index_kmer = index_kmer>>8;
    }

    return path_seq;
}

void print_uint64_kmer(uint64_t km, int km_len){
    fprintf(stdout, "uint64: %"PRIu64, km);
   fprintf(stdout, " rev: %"PRIu64, comp_rev(km, km_len));
   fprintf(stdout, " minHash: %"PRIu64, min_hash_key(km,km_len));
    unsigned char seq[33];
    seq[km_len] = '\0';
    uint64_acgt(km, seq, km_len);
    fprintf(stdout, "\tkm: %s\n", seq);

}

path_node all_nodes[Max_Path_Num * Max_Path_Len];
path_node *term_nodes[Max_Path_Num+1];

int check_kmer(uint64_t kmer){
    kmer<<48>>48;
    if(kmer > 0){
        return 1;
    }else{
        return 0;
    }
}

int find_paths(kc_c4x_t *h, unsigned char *paths[], uint64_t kmer, uint64_t mask, int cut_cov, int path_len, int k, int index_byte_num){

    all_nodes[0].kmer = kmer;
    all_nodes[0].cov = kmer_cov(min_hash_key(kmer, k),mask,h);
    int a_cov = all_nodes[0].cov/Max_Cov_Ratio, b_cov = all_nodes[0].cov * Max_Cov_Ratio ;
    if(a_cov < cut_cov + 1){a_cov = cut_cov + 1;}

    term_nodes[0] = &all_nodes[0];

    int term_node_num = 1, all_node_num = 1, p = 0;
    int seq_len = (path_len>>2) + index_byte_num;

    path_node *p_node;
    path_node next_nodes[4];
    next_nodes[0].cov = 0; next_nodes[1].cov = 0; next_nodes[2].cov = 0; next_nodes[3].cov = 0;


    path_node *p_term_nodes[Max_Path_Num+1];
    int new_term_node_num;

    while(p <= path_len && term_node_num > 0 && term_node_num <= Max_Path_Num){
        int i;
        //Copying old term nodes pointer
        for(i=0; i < Max_Path_Num; i++){
            p_term_nodes[i] = term_nodes[i]; }

        new_term_node_num = 0;

        // printf("%d\n", p);  //debuging

        int n_p = 0;
        while(n_p < term_node_num && new_term_node_num <=Max_Path_Num){

            p_node = p_term_nodes[n_p];  //This line has been modified to fix a bug

            a_cov = p_node->cov/Max_Cov_Ratio, b_cov = p_node->cov * Max_Cov_Ratio  ;
            if(a_cov < cut_cov + 1){a_cov = cut_cov + 1;}

            next_nodes[0].kmer = ((*p_node).kmer<<2) & mask;
            next_nodes[1].kmer = (((*p_node).kmer<<2) + 1ULL) & mask;
            next_nodes[2].kmer = (((*p_node).kmer<<2) + 2ULL) & mask;
            next_nodes[3].kmer = (((*p_node).kmer<<2) + 3ULL) & mask;

            next_nodes[0].cov = kmer_cov(min_hash_key(next_nodes[0].kmer, k),mask,h);
            next_nodes[1].cov = kmer_cov(min_hash_key(next_nodes[1].kmer, k),mask,h);
            next_nodes[2].cov = kmer_cov(min_hash_key(next_nodes[2].kmer, k),mask,h);
            next_nodes[3].cov = kmer_cov(min_hash_key(next_nodes[3].kmer, k),mask,h);


            int i;
            for(i=0;i<4;i++){
                if(next_nodes[i].cov >=a_cov && check_kmer(next_nodes[i].kmer) > 0 ){  
                    all_nodes[all_node_num].pre_node = p_node;
                    all_nodes[all_node_num].kmer = next_nodes[i].kmer;
                    all_nodes[all_node_num].cov = next_nodes[i].cov;

                    term_nodes[new_term_node_num] = &all_nodes[all_node_num];

                    new_term_node_num++;
                    all_node_num++;
                }
            }

            n_p++;
        }

        term_node_num = new_term_node_num;
        ++p;
    }

    int i, pass_n=0, c_index;

    for(i=0; i < term_node_num; i++){
        CALLOC(paths[i], seq_len);
        nodes_path(term_nodes[i], paths[i], path_len, k, index_byte_num);
        unsigned char *path_tmp;
        CALLOC(path_tmp, seq_len);
        int j=0;
        for(j=0; j< seq_len; j++){
            *(path_tmp + j) = one_to_two_table[*(paths[i] + j)];
        }
    }
    return term_node_num;
}


uint64_t index_encode(uint64_t index, unsigned int index_len){

    uint64_t mask = 255ULL;
    unsigned int index_byte_num = index_len>>2;
    unsigned int i;
    uint64_t m, y=0ULL;

    for(i=0; i < index_byte_num; i++){
        m = index & mask;
        m = (uint64_t)two_to_one_table[m];
        m = m<<((index_byte_num-1)*8);
        y = y>>8;
        y = y | m;
        index = index>>8;
    }
    return y;
}


uint64_t primer_uint64( char *pr, int pr_len, int index_len, int km_len){
    unsigned char *seq;
    uint64_t pr_int64;
    int i = pr_len + index_len - km_len, j=0;

    MALLOC(seq, km_len);

    while( i < pr_len ){
        *(seq + j) = *(pr + i);
        i++;
        j++;
    }

    while(j < km_len){
        *(seq + j) = 'A';
        j++;
    }

    pr_int64 = actgkmer_uint64(seq, km_len);
    free(seq);
    return pr_int64;
}

void print_kmer_seq(unsigned char *seq, int k){ //Debugging function
    int i=0;
    for(i=0;i<k;i++){
            printf("%c", *(seq + i));
    }
}

void print_cov(unsigned char* seq, int k, kc_c4x_t *h){ //Debugging function
    uint64_t mask = (1ULL<<k*2) - 1;
    int test_cov=0;
    test_cov = kmer_cov( actgkmer_hashkey(seq, k), mask, h);
    print_kmer_seq(seq, k);
    printf("\t");
    printf("cov: %d\n", test_cov);
}


int main(int argc, char *argv[])
{
    kc_c4x_t *h;
    int i, c=1, k = 31, p = KC_BITS, block_size = 10000000, n_thread = 3, min_cov = 5, index_len = 16, payload_len=128, crc_len=8, path_len;
    int max_ratio = 4;
    int mem_test = 0;

    // unsigned int min_index = 101010102, max_index=101295684, index_byte_num, max_cov=0;
    unsigned int min_index = 3111111, max_index = 3145156, index_byte_num, max_cov=0;
    uint64_t mask, pr_uint64;

    char p1[] = "CCTGCAGAGTAGCATGTC"; //p2[18] = 'CTGACACTGATGCATCCG';
    int pr_len = 18;


    ketopt_t o = KETOPT_INIT;
    while ((c = ketopt(&o, argc, argv, 1, "k:p:b:t:c:l:a:b:i:d:m:", 0)) >= 0) {
        if (c == 'k') k = atoi(o.arg);
        else if (c == 'p') p = atoi(o.arg);
        else if (c == 't') n_thread = atoi(o.arg);
        else if (c == 'c') min_cov = atoi(o.arg);
        else if (c == 'd') max_cov = atoi(o.arg);

        else if (c == 'a') min_index = atoi(o.arg);
        else if (c == 'b') max_index = atoi(o.arg);
        else if (c == 'i') index_len = atoi(o.arg);
        else if (c == 'l') payload_len = atoi(o.arg);
        else if (c == 'm') mem_test = 1;
    }

    if (argc - o.ind < 1) {
        fprintf(stderr, "\n**************************************************************************\n");
        fprintf(stderr, "**                                                                      **\n");
        fprintf(stderr, "**   DBGPS-greedy-path - De Bruijn Graph based assembler for Babel-DNA  **\n");
        fprintf(stderr, "**      Version 20230104  Author: Lifu Song lifu.song@outlook.com       **\n");
        fprintf(stderr, "**                                                                      **\n");
        fprintf(stderr, "**************************************************************************\n\n");
        fprintf(stderr, "Usage: DBGPS-greedy-path [options] <input file> \n");
        fprintf(stderr, "                       [Supporting formats: *fq, *fa, *fq.gz, *fa.gz]\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -k INT     k-mer size [%d]\n", k);
        fprintf(stderr, "  -i INT     length of index [%d] bp\n", index_len);
        fprintf(stderr, "  -l INT     data encoding length [%d] bp\n", payload_len);
        //fprintf(stderr, "  -p INT     prefix length [%d]\n", p);
        fprintf(stderr, "  -t INT     number of threads [%d]\n", n_thread);
        fprintf(stderr, "  -c INT     k-mer coverage cut-off for exclusion of noise k-mers [%d]\n", min_cov);
        fprintf(stderr, "  -d INT     Switch on k-mer coverage testing mode.  [%d]\n", max_cov);

        fprintf(stderr, "  -a INT     Initial index [%d]\n", min_index);
        fprintf(stderr, "  -b INT     End index [%d]\n", max_index);

        //fprintf(stderr, "  -r INT     Maximal ratio of k-mer coverage changes [%d]\n", max_ratio);

        fprintf(stderr, "\n");
        return 1;
    }
    if (p < KC_BITS) {
        fprintf(stderr, "ERROR: -p should be at least %d\n", KC_BITS);
        return 1;
    }


    if(max_cov < min_cov){max_cov = min_cov;}
    index_byte_num = index_len>>2;
    mask = (1ULL<<k*2) - 1;
    pr_uint64 = primer_uint64(p1, pr_len, index_len, k);

    path_len = payload_len + crc_len;

    clock_t start,end;
    double kc_dur, dec_dur;
    fprintf(stderr, "k-mers size: %d\n", k);
    fprintf(stderr, "Counting k-mers ......\n");

    start = clock();
    fprintf(stderr, "Counting file1 ......\n");
    h = count_file(argv[o.ind], k, p, block_size, n_thread);

    int f_num = argc - o.ind, c_f_n = 2;
    while(c_f_n <= f_num ) {
        fprintf(stderr, "Counting file%d ......\n", f_num);
        h = count_file2(argv[o.ind + c_f_n - 1], h, k, p, block_size, n_thread);  argv[o.ind + c_f_n - 1];
        c_f_n = c_f_n + 1;
    }

    end = clock();
    kc_dur = (double)(end - start);

    if(mem_test){
            fprintf(stderr, "Memory test mode activated! Press ENTER to continue...\n");
            getchar();}
    fprintf(stderr, "Decoding strands from index %d to ",min_index);
    fprintf(stderr, "%d ......\n",max_index);

    unsigned int id, eid;

    unsigned int cov_cut;
    int path_num;
    int dec_strand_num=0;
    for(cov_cut=min_cov;cov_cut<=max_cov;cov_cut++){
        start = clock();
        dec_strand_num = 0;

        for(id=min_index;id<=max_index;id++){
            unsigned char *path[Max_Path_Num];
            eid = index_encode(id, index_len);

            path_num = find_paths(h, path, pr_uint64 | (uint64_t)eid, mask, cov_cut, path_len, k, index_byte_num);
            if( path_num > 0){
                dec_strand_num++;
                int j=0;
                for(j=0;j<path_num;j++){
                    unsigned char DNA_seq[Max_Path_Len];
                    DNA_seq[path_len + index_byte_num*4] = '\0';
                    intseq_actg(path[j], DNA_seq, (path_len>>2) + index_byte_num);
                    printf("%d", id);

                    fprintf(stdout, "\t%s\n", DNA_seq);
                }

            }
        }

        end = clock();
        dec_dur = (double)(end - start);
        fprintf(stderr, "\nCoverage cutoff: %d ", cov_cut);
        fprintf(stderr, "\Indexes decoded: %d ", dec_strand_num);
        fprintf(stderr, "\nStrand assembly time: %f seconds\n",(dec_dur/CLOCKS_PER_SEC));

    }


    for (i = 0; i < 1<<p; ++i) kc_c4_destroy(h->h[i]);
    free(h->h); free(h);
    return 0;
}



