/* Time core_msolve without process startup or input parsing.
 *
 * Usage: msolve-inproc <input.syl> <modulus> [threads]
 *
 * Each repetition runs in a child because the installed print path omits
 * cleanup before it returns. Child exit releases that storage. A separate
 * msolve CLI call supplies the basis used for correctness checks.
 */

#define _POSIX_C_SOURCE 200809L

#include <msolve/msolve/msolve.h>

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#define MIN_REPS 3
#define MAX_REPS 5000
#define BUDGET_SECONDS 0.5
#define SINGLE_RUN_DEADLINE_SECONDS 120.0
#define TOTAL_DEADLINE_SECONDS 130.0

typedef struct {
    int32_t nvars;
    int32_t ngens;
    int32_t nterms;
    int32_t *lens;
    int32_t *exps;
    int32_t *cfs;
} syl_input;

static double now_seconds(void) {
    struct timespec timestamp;
    clock_gettime(CLOCK_MONOTONIC, &timestamp);
    return (double)timestamp.tv_sec + (double)timestamp.tv_nsec * 1e-9;
}

static void die(const char *message) {
    fprintf(stderr, "msolve-inproc: %s\n", message);
    exit(2);
}

static char *copy_string(const char *text) {
    char *copy = strdup(text);
    if (copy == NULL) {
        die("out of memory copying text");
    }
    return copy;
}

static char *read_file(const char *path) {
    FILE *file = fopen(path, "rb");
    if (file == NULL) {
        die("cannot open input file");
    }
    if (fseek(file, 0, SEEK_END) != 0) {
        die("cannot seek input file");
    }
    long size = ftell(file);
    if (size < 0 || fseek(file, 0, SEEK_SET) != 0) {
        die("cannot size input file");
    }
    char *buffer = malloc((size_t)size + 1);
    if (buffer == NULL) {
        die("out of memory reading input file");
    }
    if (fread(buffer, 1, (size_t)size, file) != (size_t)size) {
        die("short read on input file");
    }
    buffer[size] = '\0';
    fclose(file);
    return buffer;
}

static char *trim_left(char *text) {
    while (*text == ' ' || *text == '\t' || *text == '\r') {
        text++;
    }
    return text;
}

static int32_t parse_nvars(const char *text) {
    char *copy = copy_string(text);
    char *save = NULL;
    char *line = strtok_r(copy, "\n", &save);
    if (line == NULL) {
        die("empty input file");
    }
    long parsed = strtol(line, NULL, 10);
    free(copy);
    if (parsed <= 0 || parsed > INT32_MAX) {
        die("bad variable count");
    }
    return (int32_t)parsed;
}

static int32_t count_terms(char *line) {
    int32_t count = 0;
    char *save = NULL;
    for (char *term = strtok_r(line, ";", &save); term != NULL;
         term = strtok_r(NULL, ";", &save)) {
        count++;
    }
    if (count == 0) {
        die("empty generator");
    }
    return count;
}

static void count_input(const char *text, int32_t *ngens, int32_t *nterms) {
    char *copy = copy_string(text);
    char *save = NULL;
    strtok_r(copy, "\n", &save);
    *ngens = 0;
    *nterms = 0;
    for (char *line = strtok_r(NULL, "\n", &save); line != NULL;
         line = strtok_r(NULL, "\n", &save)) {
        line = trim_left(line);
        if (*line == '\0') {
            continue;
        }
        (*ngens)++;
        *nterms += count_terms(line);
    }
    free(copy);
    if (*ngens == 0) {
        die("no generators in input file");
    }
}

static syl_input allocate_input(int32_t nvars, int32_t ngens, int32_t nterms) {
    syl_input input = {
        .nvars = nvars,
        .ngens = ngens,
        .nterms = nterms,
        .lens = malloc(sizeof(int32_t) * (size_t)ngens),
        .exps = malloc(sizeof(int32_t) * (size_t)nterms * (size_t)nvars),
        .cfs = malloc(sizeof(int32_t) * (size_t)nterms),
    };
    if (input.lens == NULL || input.exps == NULL || input.cfs == NULL) {
        die("out of memory allocating msolve input arrays");
    }
    return input;
}

static void parse_term(char *term, int32_t modulus, int64_t index,
                       syl_input *input) {
    char *save = NULL;
    char *field = strtok_r(term, ",", &save);
    if (field == NULL) {
        die("empty term");
    }
    long coefficient = strtol(field, NULL, 10) % modulus;
    if (coefficient < 0) {
        coefficient += modulus;
    }
    input->cfs[index] = (int32_t)coefficient;
    for (int32_t variable = 0; variable < input->nvars; variable++) {
        field = strtok_r(NULL, ",", &save);
        if (field == NULL) {
            die("exponent vector too short");
        }
        input->exps[index * input->nvars + variable] =
            (int32_t)strtol(field, NULL, 10);
    }
    if (strtok_r(NULL, ",", &save) != NULL) {
        die("exponent vector too long");
    }
}

static int32_t parse_generator(char *line, int32_t modulus, int64_t *term_index,
                               syl_input *input) {
    int32_t count = 0;
    char *save = NULL;
    for (char *term = strtok_r(line, ";", &save); term != NULL;
         term = strtok_r(NULL, ";", &save)) {
        parse_term(term, modulus, *term_index, input);
        (*term_index)++;
        count++;
    }
    return count;
}

static void fill_input(const char *text, int32_t modulus, syl_input *input) {
    char *copy = copy_string(text);
    char *save = NULL;
    strtok_r(copy, "\n", &save);
    int32_t generator = 0;
    int64_t term = 0;
    for (char *line = strtok_r(NULL, "\n", &save); line != NULL;
         line = strtok_r(NULL, "\n", &save)) {
        line = trim_left(line);
        if (*line == '\0') {
            continue;
        }
        input->lens[generator] = parse_generator(line, modulus, &term, input);
        generator++;
    }
    free(copy);
}

static syl_input parse_syl(const char *text, int32_t modulus) {
    int32_t nvars = parse_nvars(text);
    int32_t ngens;
    int32_t nterms;
    count_input(text, &ngens, &nterms);
    syl_input input = allocate_input(nvars, ngens, nterms);
    fill_input(text, modulus, &input);
    return input;
}

static data_gens_ff_t *build_generators(syl_input input, int32_t modulus) {
    data_gens_ff_t *generators = allocate_data_gens();
    if (generators == NULL) {
        die("out of memory allocating msolve generators");
    }
    generators->nvars = input.nvars;
    generators->ngens = input.ngens;
    generators->field_char = modulus;
    generators->change_var_order = -1;
    generators->linear_form_base_coef = 0;
    generators->rand_linear = 0;
    generators->nterms = input.nterms;
    generators->lens = input.lens;
    generators->exps = input.exps;
    generators->cfs = input.cfs;
    generators->mpz_cfs = NULL;
    generators->random_linear_form =
        calloc((size_t)input.nvars, sizeof(int32_t));
    generators->vnames = malloc(sizeof(char *) * (size_t)input.nvars);
    if (generators->random_linear_form == NULL || generators->vnames == NULL) {
        die("out of memory allocating msolve metadata");
    }
    for (int32_t index = 0; index < input.nvars; index++) {
        char name[24];
        snprintf(name, sizeof name, "x%d", index + 1);
        generators->vnames[index] = copy_string(name);
    }
    generators->elim = 0;
    return generators;
}

static files_gb output_files(void) {
    files_gb files = {
        .in_file = NULL,
        .bin_file = NULL,
        .out_file = "/dev/null",
        .bin_out_file = NULL,
        .verb_file = NULL,
    };
    return files;
}

static double solve_once(data_gens_ff_t *generators, files_gb *files,
                         int32_t threads) {
    param_t *parameters = NULL;
    mpz_param_t mpz_parameters;
    mpz_param_init(mpz_parameters);
    long real_root_count = 0;
    interval *real_roots = NULL;
    real_point_t *real_points = NULL;
    double start = now_seconds();
    int status = core_msolve(
        2 /* exact sparse F4 */, 0 /* signatures */, threads, 0 /* log level */,
        17 /* initial hash size */, 0 /* pair limit */, 0 /* elimination block */,
        0 /* reset hash */, 0 /* PBM */, 1 /* reduce basis */,
        2 /* print full basis */, 0 /* truncated lifting */, 0 /* parameters */,
        2 /* genericity */, 2 /* unstable staircase */, 0 /* saturation */,
        0 /* colon */, 0 /* normal form */, 0 /* normal-form matrix */,
        0 /* basis check */, 0 /* lift matrix */, 64 /* precision */, files,
        generators, &parameters, &mpz_parameters, &real_root_count, &real_roots,
        &real_points);
    double elapsed = now_seconds() - start;
    if (status != 0) {
        fprintf(stderr, "msolve-inproc: core_msolve returned %d\n", status);
        _exit(1);
    }
    return elapsed;
}

static void child_repetition(int write_fd, data_gens_ff_t *generators,
                             files_gb *files, int32_t threads) {
    double elapsed = solve_once(generators, files, threads);
    if (write(write_fd, &elapsed, sizeof elapsed) != (ssize_t)sizeof elapsed) {
        _exit(1);
    }
    _exit(0);
}

static double run_repetition(data_gens_ff_t *generators, files_gb *files,
                             int32_t threads) {
    int descriptors[2];
    if (pipe(descriptors) != 0) {
        die("pipe failed");
    }
    pid_t child = fork();
    if (child < 0) {
        die("fork failed");
    }
    if (child == 0) {
        close(descriptors[0]);
        child_repetition(descriptors[1], generators, files, threads);
    }
    close(descriptors[1]);
    double elapsed = -1.0;
    ssize_t bytes = read(descriptors[0], &elapsed, sizeof elapsed);
    close(descriptors[0]);
    int status = 0;
    waitpid(child, &status, 0);
    if (bytes != (ssize_t)sizeof elapsed || !WIFEXITED(status) ||
        WEXITSTATUS(status) != 0) {
        die("repetition child failed");
    }
    return elapsed;
}

static int must_stop(int count, double total, double elapsed, double start) {
    if (elapsed > SINGLE_RUN_DEADLINE_SECONDS) {
        return 1;
    }
    if (now_seconds() - start > TOTAL_DEADLINE_SECONDS) {
        return 1;
    }
    return count >= MIN_REPS && total >= BUDGET_SECONDS;
}

static int run_repetitions(data_gens_ff_t *generators, files_gb *files,
                           int32_t threads) {
    double total = 0.0;
    double start = now_seconds();
    int count = 0;
    while (count < MAX_REPS) {
        double elapsed = run_repetition(generators, files, threads);
        printf("RUN %.9f\n", elapsed);
        fflush(stdout);
        total += elapsed;
        count++;
        if (must_stop(count, total, elapsed, start)) {
            break;
        }
    }
    return count;
}

static void report_usage(int count) {
    struct rusage usage;
    getrusage(RUSAGE_CHILDREN, &usage);
    printf("COUNT %d\n", count);
    printf("PEAK_RSS_KB %ld\n", usage.ru_maxrss);
}

static void parse_arguments(int argc, char **argv, int32_t *modulus,
                            int32_t *threads) {
    if (argc != 3 && argc != 4) {
        fprintf(stderr, "usage: msolve-inproc <input.syl> <modulus> [threads]\n");
        exit(2);
    }
    *modulus = (int32_t)strtol(argv[2], NULL, 10);
    *threads = argc == 4 ? (int32_t)strtol(argv[3], NULL, 10) : 1;
    if (*modulus <= 1) {
        die("modulus must be at least 2");
    }
    if (*threads < 1) {
        die("thread count must be at least 1");
    }
}

int main(int argc, char **argv) {
    int32_t modulus;
    int32_t threads;
    parse_arguments(argc, argv, &modulus, &threads);
    char *text = read_file(argv[1]);
    syl_input input = parse_syl(text, modulus);
    free(text);
    data_gens_ff_t *generators = build_generators(input, modulus);
    files_gb files = output_files();
    int count = run_repetitions(generators, &files, threads);
    report_usage(count);
    free_data_gens(generators);
    return 0;
}
