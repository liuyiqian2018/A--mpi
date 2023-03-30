/* Confidential */
/* This code is only used for CS7889d HPC@Scale at Texas State University*/

#include <cmath>
#include <vector>
#include <stdio.h>
#include <set>
#include <climits>
#include <sys/time.h>
#include <cstdlib>
#include <omp.h>
#include <mpi.h>

class BMP24 {
private:
  int wo, ho;
  int w, h;
  int* bmp;

public:
  BMP24(int xmin, int ymin, int xmax, int ymax)
  {
    if ((xmin >= xmax) || (ymin >= ymax)) exit(-2);
    wo = xmin;
    ho = ymin;
    w = xmax - xmin;
    h = ymax - ymin;
    bmp = new int[w * h];
  }

  ~BMP24()
  {
    delete [] bmp;
  }

  void dot(int x, int y, const int col)
  {
    x -= wo;
    y -= ho;
    if ((0 <= x) && (0 <= y) && (x < w) && (y < h)) {
      bmp[y * w + x] = col;
    }
  }

  void save(const char* const name)
  {
    const int pad = ((w * 3 + 3) & ~3) - (w * 3);
    FILE* f = fopen(name, "wb");
    int d;

    d = 0x4d42;  fwrite(&d, 1, 2, f);
    d = 14 + 40 + h * w * 3 + pad * h;  fwrite(&d, 1, 4, f);
    d = 0;  fwrite(&d, 1, 4, f);
    d = 14 + 40;  fwrite(&d, 1, 4, f);

    d = 40;  fwrite(&d, 1, 4, f);
    d = w;  fwrite(&d, 1, 4, f);
    d = h;  fwrite(&d, 1, 4, f);
    d = 1;  fwrite(&d, 1, 2, f);
    d = 24;  fwrite(&d, 1, 2, f);
    d = 0;  fwrite(&d, 1, 4, f);
    d = h * w * 3 + pad * h;  fwrite(&d, 1, 4, f);
    d = 0;  fwrite(&d, 1, 4, f);
    d = 0;  fwrite(&d, 1, 4, f);
    d = 0;  fwrite(&d, 1, 4, f);
    d = 0;  fwrite(&d, 1, 4, f);

    for (int y = 0; y < h; y++) {
      for (int x = 0; x < w; x++) {
        fwrite(&bmp[y * w + x], 1, 3, f);
      }
      fwrite(&d, 1, pad, f);
    }

    fclose(f);
  }
};

static int calculate_h(const int x, const int y, const int dest_x, const int dest_y, const int max_weight)
{
  return (std::abs(x - dest_x) + std::abs(y - dest_y)) * (max_weight / 2 + 1);
}

static bool on_grid(const int x, const int y, const int gsize)
{
  return (0 <= x) && (x < gsize) && (0 <= y) && (y < gsize);
}

static bool in_obst(const int idx, const int weight_n[], const int weight_e[], const int gsize)
{
  return (weight_n[idx] == INT_MAX / 2) && (weight_e[idx] == INT_MAX / 2) && (idx - 1 >= 0) && (weight_n[idx - 1] == INT_MAX / 2)
    && (idx - gsize >= 0) && (weight_e[idx - gsize] == INT_MAX / 2);
}

static bool cross(const int xa, const int ya, const int xb, const int yb, const int xc, const int yc, const int xd, const int yd)
{
  if (xa == xb || xc == xd) {
    if (xa == xb && xc == xd) {
      return (xa == xc && ((ya - yc) * (yc - yb) >= 0 || (ya - yd) * (yd - yb) >= 0));
    } else if (xa == xb) {
        return ((xc - xa) * (xa - xd) >= 0 && ((yc - ya) * (ya - yd) >= 0 || ((yc - yb) * (yb - yd) >= 0)));
    } else {
        return ((xa - xc) * (xc - xb) >= 0 && ((ya - yc) * (yc - yb) >= 0 || ((ya - yd) * (yd - yb) >= 0)));
    }
  } else {
      const float ka = (float)(yb - ya) / (xb - xa);
      const float kc = (float)(yd - yc)/ (xd - xc);
      const float da = (float)ya - ka * xa;
      const float dc = (float)yc - kc * xc;
      if (ka == kc) {
        if (da != dc) {
          return false;
        } else {
            return ((xa - xc) * (xc - xb) >= 0 || (xa - xd) * (xd - xb) >= 0);
        }
      } else {
          const float target_x = (dc - da) / (ka - kc);
          if ((xa - target_x) * (target_x - xb) >= 0 && (xc - target_x) * (target_x - xd) >= 0) {
            return true;
          } else {
            return false;
          }
      }
  }
}

static void mis(std::vector<int> fut_wl, const int num_paths, int* cur_wl, int &n_cur, const int start_idx[], const int dest_idx[], const int gsize)
{
  cur_wl[n_cur] = fut_wl[0];
  n_cur++;

  for (unsigned int i = 1; i < fut_wl.size(); i++) {
    bool flag = true;
    const int k = fut_wl[i];
    for (int j = 0; j < n_cur; j++) {
      if (cross(start_idx[k] / gsize, start_idx[k] % gsize, dest_idx[k] / gsize, dest_idx[k] % gsize, start_idx[j] / gsize, start_idx[j] % gsize, dest_idx[j] / gsize, dest_idx[j] % gsize)) {
        flag = false;
        break;
      }
    }
    if (flag) {
      cur_wl[n_cur] = k;
      n_cur++;
    }
  }
}


static void writeBMP(const int gsize, const int weight_n[], const int weight_e[], const int g[], const int start_idx[], const int dest_idx[], const std::vector<int> path[], const int num_paths, const int rank)
{
  const int red = 0x10000;
  const int green = 0x100;
  const int blue = 0x1;
  const int color_bg = 0;
  const int color_s = 255 * green;
  const int color_d = 127 * green;
  const int color_obst = 255 * red + 255 * green + 255 * blue;
  const int color_path = 255 * red;
  const int color_visited = 127 * blue;
  BMP24 bmp(0, 0, gsize, gsize);

  for (int i = 0; i < gsize * gsize; i++) {
    const int x = i / gsize;
    const int y = i % gsize;
    int col = color_bg;
    if (in_obst(i, weight_n, weight_e, gsize)) {
      col = color_obst;
    } else if (g[i] < INT_MAX / 2) {
        col = color_visited;
    }
    bmp.dot(x, y, col);
  }

  for (int i = 0; i < num_paths; i++) {
    std::vector<int> tpath = path[i];
    for (int j = 0; j < (int)tpath.size(); j++) {
      const int x = tpath[j] / gsize;
      const int y = tpath[j] % gsize;
      bmp.dot(x, y, color_path);
    }
  }

  for (int i = 0; i < num_paths; i++) {
    const int start_x = start_idx[i] / gsize;
    const int start_y = start_idx[i] % gsize;
    const int dest_x = dest_idx[i] / gsize;
    const int dest_y = dest_idx[i] % gsize;
    bmp.dot(start_x, start_y, color_s);
    bmp.dot(dest_x, dest_y, color_d);
  }
  char name[256];
  sprintf(name, "picture%d.bmp", rank);

  bmp.save(name);
}

static std::vector<int> trace_path(const int gsize, const int weight_n[], const int weight_e[], const int g[], const int start_idx, const int dest_idx)
{
  const int offset[4] = {-1, -gsize, 1, gsize};
  std::vector<int> path;
  int current_idx = dest_idx;
  int parent_idx = current_idx + offset[g[current_idx] % 4];
  int previous_idx = current_idx;

  path.push_back(current_idx);
  int len = 2;
  int turns = 0;

  while ((parent_idx != start_idx) && (previous_idx != parent_idx) && (g[current_idx] < INT_MAX / 2)) {
    len++;
    if ((len != 3) && (g[parent_idx] % 4 != g[current_idx] % 4)) {
      turns++;
    }
    path.push_back(parent_idx);
    previous_idx = current_idx;
    current_idx = parent_idx;
    parent_idx = current_idx + offset[g[current_idx] % 4];
  }
  path.push_back(parent_idx);

  if (parent_idx != start_idx) {
    path.clear();
  } else {
    /*printf("path length: %d\n", len);*/
    /*printf("turns: %d\n", turns);*/
    /*printf("g_start: %d\n", g[start_idx] / 4);*/
    /*printf("g_dest: %d\n\n", g[dest_idx] / 4);*/
  }
  return path;
}

static void find_path(const int gsize, const int weight_n[], const int weight_e[], int g[], const int start_idx, const int dest_idx, const int penalty, const int baseval, const int max_weight)
{
  /*printf("baseval: %d, g_start: %d\n", baseval, g[start_idx] / 4);*/
  if (baseval < 0) return;

  int t;
  do {
    #pragma omp atomic read
    t = g[start_idx];
    if (t / 4 < baseval) return;
  } while (__sync_val_compare_and_swap(&g[start_idx], t, baseval * 4) != t);

  std::set< std::pair<int, int> > open_set;

  const int dest_x = dest_idx / gsize;
  const int dest_y = dest_idx % gsize;
  const int temp_h = calculate_h(start_idx / gsize, start_idx % gsize, dest_x, dest_y, max_weight);
  const int temp_g = baseval * 4;
  open_set.insert(std::make_pair(temp_h + temp_g, start_idx));

  int iter = 0;
  while (!open_set.empty()) {
    iter++;
    const int current_idx = (*open_set.begin()).second;
    open_set.erase(open_set.begin());

    if (current_idx == dest_idx) {
      /*printf("iterations: %d\n", iter);*/
      return;
    }

    const int current_x = current_idx / gsize;
    const int current_y = current_idx % gsize;
    int weight_arr[4] = {weight_n[current_idx], weight_e[current_idx], INT_MAX / 2, INT_MAX / 2};
    if (current_y > 0) {
      weight_arr[2] = weight_n[current_x * gsize + current_y - 1];
    }
    if (current_x > 0) {
      weight_arr[3] = weight_e[(current_x - 1) * gsize + current_y];
    }

    const int x_offset[4] = {0, 1, 0, -1};
    const int y_offset[4] = {1, 0, -1, 0};

    for (int i = 0; i < 4; i++) {
      const int nbr_x = current_x + x_offset[i];
      const int nbr_y = current_y + y_offset[i];
      const int nbr_idx = nbr_x * gsize + nbr_y;

      if (on_grid(nbr_x, nbr_y, gsize) && !in_obst(nbr_idx, weight_n, weight_e, gsize)) {
        int temp;
        #pragma omp atomic read
        temp = g[current_idx];
        const int gc = temp;
        #pragma omp atomic read
        temp = g[nbr_idx];
        const int gn = temp;
        if (gc / 4 < baseval) return;
        if (gn / 4 < baseval) return;

        const int temp_h = calculate_h(nbr_x, nbr_y, dest_x, dest_y, max_weight);
        int temp_g = gc / 4 + weight_arr[i];

        if (gc % 4 != i) {
          temp_g += penalty;
        }

        if (temp_g < gn / 4) {

          do {
          #pragma omp atomic read
          t = g[nbr_idx];
          if (t / 4 < baseval) return;
          } while (__sync_val_compare_and_swap(&g[nbr_idx], t, temp_g * 4 + i) != t);
          open_set.insert(std::make_pair(temp_g + temp_h, nbr_idx));
        }
      }
    }
  }
}

int main(int argc, char *argv[])
{
  if (argc != 5) {
    fprintf(stderr, "USAGE: %s grid_size num_paths num_obstacles penalty\n", argv[0]);
    exit(-1);
  }

  int nproc, rank;
  int num_elements_per_proc = atoi(argv[1]);

  MPI_Init (&argc, &argv); /* Initialize MPI */
  MPI_Comm_size(MPI_COMM_WORLD,&nproc); /* Get Comm Size - how many ppl working*/
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); /* Get rank -- what is my position */

  const int gsize = atoi(argv[1]);
  const int num_paths = atoi(argv[2]) / nproc;
  const int num_obstacles = atoi(argv[3]);
  const int penalty = atoi(argv[4]);
  const int max_weight = 10;

  srand(rank);
  int* const weight_n = new int [gsize * gsize];
  int* const weight_e = new int [gsize * gsize];
  int* const g = new int [gsize * gsize];
  for (int i = 0; i < gsize; i++) {
    for (int j = 0; j < gsize; j++) {
      const int k = i * gsize + j;
      weight_n[k] = rand() % max_weight + 1;
      weight_e[k] = rand() % max_weight + 1;
      g[k] = INT_MAX / 2;
    }
  }

  const int osize = gsize / 16;
  for (int i = 0; i < num_obstacles; i++) {
    const int obst_x = rand() % (gsize - 1 - osize) + 1;
    const int obst_y = rand() % (gsize - 1 - osize) + 1;
    for (int i = obst_x; i < obst_x + osize; i++) {
      for (int j = obst_y; j < obst_y + osize; j++) {
        const int k = i * gsize + j;
        weight_n[k] = INT_MAX / 2;
        weight_e[k] = INT_MAX / 2;
        weight_n[k - 1] = INT_MAX / 2;
        weight_e[k - gsize] = INT_MAX / 2;
      }
    }
  }

  int* start_idx = new int [num_paths];
  int* dest_idx = new int [num_paths];
  for (int i = 0; i < num_paths; i++) {
    do {
      const int start_x = rand() % gsize;
      const int start_y = rand() % gsize;
      const int dest_x = rand() % gsize;
      const int dest_y = rand() % gsize;
      start_idx[i] = start_x * gsize + start_y;
      dest_idx[i] = dest_x * gsize + dest_y;
    } while ((start_idx[i] == dest_idx[i]) || in_obst(start_idx[i], weight_n, weight_e, gsize) || in_obst(dest_idx[i], weight_n, weight_e, gsize));
  }

  std::vector<int>* path = new std::vector<int> [num_paths];
  const int gbase = gsize * max_weight * 4;
  int base = num_paths * num_paths;
  const int size = sqrt(gsize / 10);
  bool* found = new bool [num_paths];
  std::fill(found, found + num_paths, false);

  int* cur_wl = new int [num_paths];
  int n_cur = 0;
  std::vector<int> fut_wl;
  for (int i = 0; i < num_paths; i++) {
    fut_wl.push_back(i);
  }

  timeval start, end;
  gettimeofday(&start, NULL);

  while (!fut_wl.empty() && base >= 0) {
    mis(fut_wl, num_paths, cur_wl, n_cur, start_idx, dest_idx, gsize);

    #pragma omp parallel for default(none) shared(cur_wl, n_cur, fut_wl, found, start_idx, dest_idx, base, path, gbase, penalty, g, weight_e, weight_n, gsize) schedule(dynamic)
    for (int i = 0; i < n_cur; i++) {
      const int k = cur_wl[i];
      find_path(gsize, weight_n, weight_e, g, start_idx[k], dest_idx[k], penalty, (base + k) * gbase, max_weight);
      path[k] = trace_path(gsize, weight_n, weight_e, g, start_idx[k], dest_idx[k]);
      if (path[k].size() > 0) {
        found[k] = true;
      }
    }

    const int p = n_cur;
    const int c = fut_wl.size();
    n_cur = 0;
    base -= num_paths;
    fut_wl.clear();
    for (int i = 0; i < num_paths; i++) {
      if (!found[i]) {
        fut_wl.push_back(i);
      }
    }
    const int n = c - fut_wl.size();
    /* printf("number of paths found in this iteraion: %d\n\n", n);
    printf("number of paths not found in this iteraion: %d\n\n", (p - n)); */

  }

  MPI_Barrier(MPI_COMM_WORLD);
  gettimeofday(&end, NULL);
  if (rank == 0) {
    const double runtime = end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("compute time: %.4f s\n", runtime);
  }

  int count = 0;
  for (int i = 0; i < num_paths; i++) {
    if (found[i]) count++;
  }
  printf("%d found sucessful paths: %d\n", rank, count);

  writeBMP(gsize, weight_n, weight_e, g, start_idx, dest_idx, path, num_paths, rank);

  delete [] weight_n;
  delete [] weight_e;
  delete [] g;
  delete [] start_idx;
  delete [] dest_idx;
  delete [] path;
  MPI_Finalize(); /* Finalize */
  return 0;
}
