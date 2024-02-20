#define _CRT_SECURE_NO_WARNINGS
#include<cstring>
#include<map>
#include<vector>
#include<cmath>
#include<algorithm>
using namespace std;

struct iontype {
   char s[10];//이온식
   int n;//개수
   int ind;//전하량
};
struct subtype {
   char s[10];//화학식
   int n;//개수
   iontype ion[3];//ion[1] : 양이온, ion[2] : 음이온
};
struct n4 {
   int coef[5];
};

char input[3][50];//입력받는 식
char s[10];
subtype sub[3][3];
map<int, n4> m;
vector <int> v;
char ends[6][4];//각 화학식 뒤에 출력할 문자열(+, -> 등)
void PolyatomicIon(int rep, int index, int check, int i, int j, int len) { //다원자 이온 정리, 개수 구하기
   if (check == 1) {
      sub[rep][i].ion[j].n = s[index + len + 1] - '0';
      sub[rep][i].ion[j].s[index + len + 1] = '\0';
   }
   else sub[rep][i].ion[j].n = 1;
   return;
}

void FindIonNum(int rep, int index, int i, int j) { //일반 이온의 개수 구하기
   for (; s[index] != '\0'; index++);
   if ('2' <= s[index - 1] && s[index - 1] <= '9') {
      sub[rep][i].ion[j].n = s[index - 1] - '0';
      sub[rep][i].ion[j].s[index - 1] = '\0';
   }
   else sub[rep][i].ion[j].n = 1;
   return;
}
void DeleteBracket(int rep, int i, int j) { //괄호 제거
   int index;
   for (index = 1; s[index] != ')'; index++) {
      sub[rep][i].ion[j].s[index - 1] = s[index];
      sub[rep][i].ion[j].s[index] = '\0';
   }
   sub[rep][i].ion[j].s[index] = '\0';
   return;
}
void Seperate(int rep, int index, int i, int j) { // 양이온, 음이온 분리
   for (int k = index; k < j; k++) {
      sub[rep][i].ion[1].s[k - index] = s[k];
   }
   for (int k = 0; s[j + k] != '\0'; k++) {
      sub[rep][i].ion[2].s[k] = s[j + k];
   }
   return;
}

void ResetS(int rep, int i, int k) { //s 재설정
   for (int j = 0; j < 10; j++) {
      s[j] = '\0';
   }
   for (int j = 0; sub[rep][i].s[j] != '\0'; j++) {
      s[j] = sub[rep][i].ion[k].s[j];//s 재설정
   }
   return;
}
void FindIndex(int rep, int i, int j) {
   for (int k = 0; sub[rep][i].ion[j].s[k] != '\0'; k++) {
      s[k] = sub[rep][i].ion[j].s[k];
   }
   if (strncmp(s, "Li", 2) == 0) sub[rep][i].ion[j].ind = 3;
   else if (strncmp(s, "Na", 2) == 0) sub[rep][i].ion[j].ind = 4;
   else if (strncmp(s, "Ag", 2) == 0) sub[rep][i].ion[j].ind = 6;
   else if (strncmp(s, "NH4", 3) == 0) sub[rep][i].ion[j].ind = 7;
   else if (strncmp(s, "Be", 2) == 0) sub[rep][i].ion[j].ind = 8;
   else if (strncmp(s, "Mg", 2) == 0) sub[rep][i].ion[j].ind = 9;
   else if (strncmp(s, "Ca", 2) == 0) sub[rep][i].ion[j].ind = 10;
   else if (strncmp(s, "Ba", 2) == 0) sub[rep][i].ion[j].ind = 11;
   else if (strncmp(s, "Pb", 2) == 0) sub[rep][i].ion[j].ind = 12;
   else if (strncmp(s, "Cu", 2) == 0) sub[rep][i].ion[j].ind = 13;
   else if (strncmp(s, "Zn", 2) == 0) sub[rep][i].ion[j].ind = 14;
   else if (strncmp(s, "Al", 2) == 0) sub[rep][i].ion[j].ind = 15;
   else if (strncmp(s, "Fe", 2) == 0) sub[rep][i].ion[j].ind = 16;
   else if (strncmp(s, "Cl", 2) == 0) sub[rep][i].ion[j].ind = 18;
   else if (strncmp(s, "Br", 2) == 0) sub[rep][i].ion[j].ind = 19;
   else if (strncmp(s, "OH", 2) == 0) sub[rep][i].ion[j].ind = 21;
   else if (strncmp(s, "NO3", 3) == 0) sub[rep][i].ion[j].ind = 22;
   else if (strncmp(s, "CO3", 2) == 0) sub[rep][i].ion[j].ind = 25;
   else if (strncmp(s, "SO4", 2) == 0) sub[rep][i].ion[j].ind = 26;
   else if (strncmp(s, "C", 1) == 0) sub[rep][i].ion[j].ind = 1;
   else if (strncmp(s, "H", 1) == 0) sub[rep][i].ion[j].ind = 2;
   else if (strncmp(s, "K", 1) == 0) sub[rep][i].ion[j].ind = 5;
   else if (strncmp(s, "F", 1) == 0) sub[rep][i].ion[j].ind = 17;
   else if (strncmp(s, "I", 1) == 0) sub[rep][i].ion[j].ind = 20;
   else if (strncmp(s, "O", 1) == 0) sub[rep][i].ion[j].ind = 23;
   else if (strncmp(s, "S", 1) == 0) sub[rep][i].ion[j].ind = 24;
   return;
}
int gcd(int a, int b)
{
   while (b != 0)
   {
      int temp = a % b;
      a = b;
      b = temp;
   }
   return a;
}
int subisone[3] = { 0, 1, 1 };//물질이 하나인지
int main() {
   printf("화학식을 입력하세요\n");
   scanf("%s", input[0]);
   //좌변 우변 ”->“ 좌우로 분리
   for (int i = 0; input[0][i] != '\0'; i++) {
      if (input[0][i] == '-') {
         for (int j = 0; j < i; j++) {
            input[1][j] = input[0][j];
         }
         for (int j = 0; input[0][i + 2 + j]; j++) {
            input[2][j] = input[0][i + 2 + j];
         }
         break;
      }
   }
   for (int rep = 1; rep <= 2; rep++)
   {
      //각 변의 식 ’+‘ 좌우로 분리
      for (int i = 0; input[rep][i] != '\0'; i++) {
         if (input[rep][i] == '+') {
            for (int j = 0; j < i; j++) {
               sub[rep][1].s[j] = input[rep][j];
            }
            for (int j = 0; input[rep][i + 1 + j] != '\0'; j++) {
               sub[rep][2].s[j] = input[rep][i + 1 + j];
            }
            subisone[rep] = 0;
            break;
         }
      }
      if (subisone[rep] == 1)
      {   
         for (int j = 0; input[rep][j] != '\0'; j++) {
            sub[rep][1].s[j] = input[rep][j];
         }
      }
      for (int i = 1; i <= 2; i++) {
         int check = 0;//괄호가있는지없는지
         for (int j = 0; j < 10; j++) {
            s[j] = '\0';
         }
         for (int j = 0; sub[rep][i].s[j] != '\0'; j++) {
            s[j] = sub[rep][i].s[j];
         }
         int index = 0;

         //양이온, 음이온 분리
         int ifionisone = 1;
         for (int j = index; s[j] != '\0'; j++) {
            if (j != index) {
               if (s[j] == '(') {
                  ifionisone = 0;
                  Seperate(rep, index, i, j);
                  break;
               }
               else if (('A' <= s[j] && s[j] <= 'Z') && !('A' <= s[j - 1] && s[j - 1] <= 'Z')) {
                  ifionisone = 0;
                  Seperate(rep, index, i, j);
                  break;
               }
            }
         }
         if (ifionisone == 1) { //이온이 2개가 아닐 경우 
            for (int j = index; s[j] != '\0'; j++) {
               sub[rep][i].ion[1].s[j] = s[j];
            }
         }

         index = 0;//양이온 분석 시작
         ResetS(rep, i, 1);
         //다원자이온 체크
         check = 0;
         if (s[0] == '(') {
            index += 1;
            check = 1;
         }

         if (s[index] == 'N' && s[index + 1] == 'H' && s[index + 2] == '4') {//다원자이온(NH4)
            PolyatomicIon(rep, index, check, i, 1, 3);
         }
         else { // 다원자이온이 아니라면
            FindIonNum(rep, index, i, 1);
         }

         //괄호 제거
         if (check == 1) {
            DeleteBracket(rep, i, 1);
         }

         //번호 붙이기
         FindIndex(rep, i, 1);

         if (ifionisone == 1) { // 이온 2개가 아닌 경우
            continue;
         }
         //———————————————————————————————————————
         index = 0;//음이온 분석 시작
         ResetS(rep, i, 2);

         //다원자이온 체크
         check = 0;
         if (s[0] == '(') {
            index += 1;
            check = 1;
         }

         if (s[index] == 'N' && s[index + 1] == 'O' && s[index + 2] == '3') { //다원자이온(NO3)
            PolyatomicIon(rep, index, check, i, 2, 3);
         }
         else if (s[index] == 'C' && s[index + 1] == 'O' && s[index + 2] == '3') { //다원자이온(CO3)
            PolyatomicIon(rep, index, check, i, 2, 3);
         }
         else if (s[index] == 'S' && s[index + 1] == 'O' && s[index + 2] == '4') { //다원자이온(SO4)
            PolyatomicIon(rep, index, check, i, 2, 3);
         }
         else if (s[index] == 'O' && s[index + 1] == 'H') { //다원자이온(OH)
            PolyatomicIon(rep, index, check, i, 2, 2);
         }
         else { //다원자이온이 아니라면
            FindIonNum(rep, index, i, 2);
         }

         //괄호 제거
         if (check == 1) {
            DeleteBracket(rep, i, 2);
         }

         //번호 붙이기
         FindIndex(rep, i, 2);
      }
   }
   int n[5] = {0, 0 ,0, 0, 0};
   int subisexist[5] = { 0, 1, 1, 1, 1 };
   if (subisone[1] == 1) subisexist[2] = 0;
   if (subisone[2] == 1) subisexist[4] = 0;
   n[1] = 2520;//1~10의 최소공배수, 분수 저장이 불가능하기에 정수비로 표현
   //연립방정식 작성
   for (int i = 1; i <= 2; i++)
   {
      for (int j = 1; j <= 2; j++)
      {
         for (int k = 1; k <= 2; k++)
         {
            if (sub[i][j].ion[k].n == 0) continue;
            m[sub[i][j].ion[k].ind].coef[(i - 1) * 2 + j] = sub[i][j].ion[k].n;
            if (i == 1) v.push_back(sub[i][j].ion[k].ind);
         }
      }
   }

   //연립방정식 풀이
   int s = v.size();
   for (int rep = 1; rep <= 4; rep++)
   {
      for (int i = 0; i < s; i++)
      {
         int check = 0;//각 등식에서 미지수 한개만 모르는지 체크 : m[i].n[j] != 0 && n[j] == 0
         int xind = 0;
         for (int j = 1; j <= 4; j++)
         {
            if (subisexist[j] == 1 && m[v[i]].coef[j] != 0 && n[j] == 0)
            {
               check++;
               xind = j;
            }
         }
         if (check == 1)
         {
            n[xind] = abs((m[v[i]].coef[1] * n[1] + m[v[i]].coef[2] * n[2]) - (m[v[i]].coef[3] * n[3] + m[v[i]].coef[4] * n[4])) / m[v[i]].coef[xind];
         }
      }
   }

   

   //정수비 간단하게 정리
   int G = gcd(gcd(n[1], n[2]), gcd(n[3], n[4]));
   for (int i = 1; i <= 2; i++)
   {
      for (int j = 1; j <= 2; j++)
      {
         sub[i][j].n = n[(i - 1) * 2 + j] / G;
      }
   }
   if(subisexist[2] == 1) ends[1][0] = '+'; 
   ends[2][0] = '-'; ends[2][1] = '>'; 
   if (subisexist[4] == 1) ends[3][0] = '+';  //각 화학식 뒤에 출력할 문자열 지정(+, ->)
   for (int i = 1; i <= 2; i++)
   {
      for (int j = 1; j <= 2; j++)
      {
         if (sub[i][j].n > 1)
         {
            printf("%d", sub[i][j].n);
         }
         printf("%s%s", sub[i][j].s, ends[(i - 1) * 2 + j]);
      }
   }
}
