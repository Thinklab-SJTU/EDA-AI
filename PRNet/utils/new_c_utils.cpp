#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <set>

using namespace std;

// To compile:
// gcc utils.c -shared -fPIC -o utils.so
// g++ -std=c++11 c_utils.cpp -shared -fPIC -o c_utils.so

extern "C"{
 
    int process(double * grid, double * A, int m, int n, int * num, bool * over);
 
}


bool dfs(double * grid, int r, int c, int num, int m, int n, set<int>& pin, set<int>& record, bool & over) {
    bool flag = false;
    bool temp_over = false;
    grid[r * n + c] = num;
    if (pin.count(r * n + c)) flag = true;
    record.erase(r * n + c);
    if (record.empty()) {
        over = true;
        return true;
    }

    int a[4] = {-1, 1, 0, 0};
    int b[4] = {0, 0, -1, 1};

    for (int i = 0; i < 4; i++) {
        int x = r + a[i], y = c + b[i];
        if (0 <= x && x < m && 0 <= y && y < n && grid[x * n + y] == 1) {
            bool temp_flag = dfs(grid, x, y, num, m, n, pin, record, temp_over);
            flag = flag || temp_flag;
            over = over || temp_over;
            if (over) return true;
        }
    }

    if (!pin.count(r * n + c) && !flag) {
        grid[r * n + c] = 0;
        return false;
    }

    return true;
}

int process(double * grid, double * A, int m, int n, int * num, bool * over) {

    for (int i = 0; i < m * n; i++) grid[i] = (grid[i] > 0) ? 1 : 0;

    int index[m * n];
    int cnt = 0;
    set<int> pin;
    int hpwl = 0;
    int max_x = 0, max_y = 0;
    int min_x = m + n, min_y = m + n;

    for (int i = 0; i < m * n; i++) {
        if (A[i] > 0) {
            max_x = max(max_x, i / n);
            max_y = max(max_y, i % n);
            min_x = min(min_x, i / n);
            min_y = min(min_y, i % n);
            index[cnt] = i;
            cnt++;
            pin.insert(i);
        }
    }

    set<int> record(pin);
    hpwl = max_x - min_x + max_y - min_y - 1;

    for (int i = 0; i < cnt; i++) {
        int x = index[i] / n;
        int y = index[i] % n;
        if (grid[x * n + y] == 1) {
            bool flag;
            flag = dfs(grid, x, y, *num, m, n, pin, record, *over);
            if (*over) break;
            *num += 1;
        }
    }

    for (int i = 0; i < m * n; i++) grid[i] = (grid[i] <= 1) ? 0 : 1;

    return hpwl;
}