/*
 * MacQC C API.
 *
 *   macqc_run(input_json)  runs HF (+ FCI) and returns the result as a JSON string.
 *                          It never throws and never returns NULL except when memory
 *                          is exhausted. Errors are reported in the "error" field.
 *   macqc_free(p)          releases a string returned by macqc_run.
 */
#ifndef MACQC_API_H
#define MACQC_API_H

#ifdef __cplusplus
extern "C" {
#endif

char* macqc_run(const char* input_json);
void macqc_free(char* p);

#ifdef __cplusplus
}
#endif

#endif /* MACQC_API_H */
