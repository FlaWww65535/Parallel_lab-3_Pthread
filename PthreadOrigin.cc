#include<iostream>
#include<cmath>
#include<pthread.h>
#include<semaphore.h>
#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>

using namespace std;


#define NUM_THREADS 7

const int MAXN=2048;
int N;
const double eps=1e-3;

float a[MAXN][MAXN];
float b[MAXN][MAXN];

char cases[6][20]={"input64.dat","input128.dat","input256.dat","input512.dat","input1024.dat","input2048.dat"};
int range[6]={64,128,256,512,1024,2048};

class CTimer
{
public:
	CTimer(void);
	~CTimer(void);

	void time_in();
	long long time_out();

private:
	struct timeval start;
	struct timeval end;
};

CTimer::CTimer(void)
{
}


CTimer::~CTimer(void)
{
}

void CTimer::time_in()
{
	gettimeofday(&start,NULL);
}

long long CTimer::time_out()
{
	gettimeofday(&end,NULL);
	long long timer=1000000 * (end.tv_sec-start.tv_sec)+ end.tv_usec-start.tv_usec;
	return timer;
}

void read(){
	for(int i=0;i<N;i++)
		for(int j=0;j<N;j++){
			cin>>a[i][j];
		}
	for(int i=0;i<N;i++)
		for(int j=0;j<N;j++){
			cin>>b[i][j];
		}
}


class Trivial{
public:
	float A[MAXN][MAXN];
	float res[MAXN][MAXN];
	void read(){
		for(int i=0;i<N;i++)
			for(int j=0;j<N;j++){
				A[i][j]=a[i][j];
			}
		for(int i=0;i<N;i++)
			for(int j=0;j<N;j++){
				res[i][j]=b[i][j];
			}
	}
	void calculate(){
		CTimer ct;
		ct.time_in();
		for(int k=0;k<N;k++){
			for(int j=k+1;j<N;j++){
				A[k][j]=A[k][j]/A[k][k];
			}
			A[k][k]=1.0f;
			for(int i=k+1;i<N;i++){
				for(int j=k+1;j<N;j++){
					A[i][j]=A[i][j]-A[i][k]*A[k][j];
				}
				A[i][k]=0;
			}
		}
		printf("trivial algo costs:%lldus",ct.time_out());
	}
	void check(){
		bool flag=0;
		for(int i=0;i<N;i++)
			for(int j=0;j<N;j++){
				if(abs(A[i][j]-res[i][j])>eps)flag=1;
			}
		if(flag){
			cout<<"wrong"<<endl;
		}
		else{
			cout<<"correct"<<endl;
		}
	}

}trivial_algo;

typedef struct {
	int t_id; //线程 id
}threadParam_t;

 //信号量定义
sem_t sem_main;
sem_t sem_workerstart[NUM_THREADS]; // 每个线程有自己专属的信号量
sem_t sem_workerend[NUM_THREADS];

static void *threadFunc(void *param) {
	threadParam_t *p = (threadParam_t*)param;
	int t_id = p -> t_id;
	for(int k=0;k<N;++k){
		sem_wait(&sem_workerstart[t_id]);
		for(int i=k+1+t_id; i < N; i += NUM_THREADS){
			for(int j=k+1;j<N;++j){
				a[i][j] = a[i][j]-a[i][k]*a[k][j];
			}
			a[i][k]=0.0;
		}
		sem_post(&sem_main);
		sem_wait(&sem_workerend[t_id]);
	}
	pthread_exit(NULL);
}
void destroySem(){
	for(int t_id = 0; t_id < NUM_THREADS; t_id++){
		sem_destroy(&sem_workerstart[t_id]);
		sem_destroy(&sem_workerend[t_id]);
	}
	sem_destroy(&sem_main);
}
class Pthread{
public:
	float res[MAXN][MAXN];
	void read(){

		for(int i=0;i<N;i++)
			for(int j=0;j<N;j++){
				res[i][j]=b[i][j];
			}
	}

	void calculate(){
		CTimer ct;
		ct.time_in();
		sem_init(&sem_main, 0, 0);
		for(int i = 0; i < NUM_THREADS; ++i){
			sem_init(&sem_workerstart[i], 0, 0);
			sem_init(&sem_workerend[i], 0, 0);
		}
		//创建线程
		pthread_t handles[NUM_THREADS];// 创建对应的 Handle
		threadParam_t param[NUM_THREADS];// 创建对应的线程数据结构
		for(int t_id = 0; t_id < NUM_THREADS; t_id++){
			param[t_id].t_id = t_id;
			pthread_create(handles+t_id,NULL,threadFunc,param+t_id);
		}
		for(int k = 0; k < N; ++k){
			//主线程做除法操作
			for(int j = k+1; j < N; j++){
				a[k][j] = a[k][j] / a[k][k];
			}
			a[k][k] = 1.0;
			//开始唤醒工作线程
			for (int t_id = 0; t_id < NUM_THREADS; ++t_id){
				sem_post(&sem_workerstart[t_id]);
			}
			//主线程睡眠（等待所有的工作线程完成此轮消去任务）
			for (int t_id = 0; t_id < NUM_THREADS; ++t_id){
				sem_wait(&sem_main);
			}
			// 主线程再次唤醒工作线程进入下一轮次的消去任务
			for (int t_id = 0; t_id < NUM_THREADS; ++t_id){
				sem_post(&sem_workerend[t_id]);
			}
		}

		for(int t_id = 0; t_id < NUM_THREADS; t_id++){
			pthread_join(handles[t_id],NULL);
		}
		//销毁所有信号量
		destroySem();
		printf("Pthread algo costs:%lldus",ct.time_out());
	}
	void check(){
		bool flag=0;
		for(int i=0;i<N;i++)
			for(int j=0;j<N;j++){
				if(abs(a[i][j]-res[i][j])>eps){if(flag==0)cout<<a[i][j]<<" "<<res[i][j];flag=1;}
			}
		if(flag){
			cout<<"wrong"<<endl;
		}
		else{
			cout<<"correct"<<endl;
		}
	}
}Pthread_algo;



int main() {
	for(int i=0;i<6;i++){
		N=range[i];
		freopen(cases[i],"r",stdin);
		cout<<cases[i]<<" condition:"<<endl;
		read();
		trivial_algo.read();
		trivial_algo.calculate();
		trivial_algo.check();
		Pthread_algo.read();
		Pthread_algo.calculate();
		Pthread_algo.check();
	}
	return 0;
}
