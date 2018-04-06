#ifndef __PROGTEST__

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <deque>
#include <queue>
#include <stack>
#include <algorithm>
#include <pthread.h>
#include <semaphore.h>
#include <cstdint>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <mutex>
#include <memory>
#include <condition_variable>
#include <atomic>
using namespace std;


class CFITCoin;
class CCVUTCoin;
class CCustomer;

typedef struct shared_ptr<CFITCoin>                        AFITCoin;
typedef struct shared_ptr<CCVUTCoin>                       ACVUTCoin;
typedef struct shared_ptr<CCustomer>                       ACustomer;
//=================================================================================================
class CFITCoin
{
  public:
                             CFITCoin                      ( const vector<uint32_t> & vectors,
                                                             int               distMax )
                             : m_Vectors ( vectors ),
                               m_DistMax ( distMax ),
                               m_Count ( 0 )
    {
    }
    virtual                  ~CFITCoin                     ( void ) = default;
    vector<uint32_t>         m_Vectors;
    int                      m_DistMax;
    uint64_t                 m_Count;
};
//=================================================================================================
class CCVUTCoin       
{
  public:
                             CCVUTCoin                     ( const vector<uint8_t> & data,
                                                             int               distMin,
                                                             int               distMax )
                             : m_Data ( data ),
                               m_DistMin ( distMin ),
                               m_DistMax ( distMax ),
                               m_Count ( 0 )
    {
    }
    virtual                  ~CCVUTCoin                    ( void ) = default;
    vector<uint8_t>          m_Data;
    int                      m_DistMin;
    int                      m_DistMax;
    uint64_t                 m_Count;
};
//=================================================================================================
class CCustomer
{
  public:
    virtual                  ~CCustomer                    ( void ) = default;
    virtual AFITCoin         FITCoinGen                    ( void ) = 0;
    virtual ACVUTCoin        CVUTCoinGen                   ( void ) = 0;
  
    virtual void             FITCoinAccept                 ( AFITCoin          x ) = 0;
    virtual void             CVUTCoinAccept                ( ACVUTCoin         x ) = 0;
};
//=================================================================================================
#endif /* __PROGTEST__ */

int BUFFER_SIZE_CRIG = 40;
int BUFFER_SIZE_CUSTOMER = 20;

class CRig
{
  private:
    class Customer;
  public:
    static void Solve(ACVUTCoin x){
        string str = "";
        vector<uint8_t>& vec = x->m_Data;
        ConvertToString(vec, str);
        // pouze pro debugovani
        //cout << str << endl;
        x->m_Count = GeneratePrefixSufix(str, x->m_DistMin, x->m_DistMax);
    }
    static void Solve(AFITCoin x){
        if (x->m_Vectors.size() == 0){
            return;
        }
        vector<uint32_t> vec;
        uint32_t fix = 0;
        // SimplifiedVector
        fix = SimplifiedVector(x->m_Vectors, vec);
        x->m_Count = FindDist(vec, x->m_DistMax, fix);
    }
    CRig(void): stop_flag(false) {
        sem_init(&sem_free, 0, BUFFER_SIZE_CRIG);
        sem_init(&sem_full, 0, 0);
    }
    ~CRig(void){
        sem_destroy(&sem_free);
        sem_destroy(&sem_full);

        for (Customer* c : customers){
            delete c;
        }
    }
    void Start(int thrCnt){
        this->thrCnt = thrCnt;
        // spusti pracovnich thrCnt vlaken a bude odebirat a obsluhovat zakazniky z fronty nebo vectoru
        //cout << "Start running" << endl;
        for ( int i = 0; i < thrCnt; ++i){
            workers.emplace_back(&CRig::WorkingThread, this);
        }
    }
    void Stop(void){
        // joinovat vlakna
        for (auto& t : producents){
            t.join();
        }
        {
            unique_lock<mutex> ul(mtx_stop);
            stop_flag = true;
        }
        // podle poctu vlaken volam pocet probuzeni
        for (int i = 0; i < thrCnt; ++i){
            sem_post(&sem_full);
        }
        
        for (auto& t : workers){
            t.join();
        }
        // mozna jeste sem pridat postovani na sem_fill !!!
        for (auto& c : customers){
            sem_post(&c->sem_fill);
        }

        for (auto& t : consumers){
            t.join();
        }

        /*cout << "buffer size: " << buffer.size() << endl;
        for (const auto& c : customers)
            cout << "   done size: " << c->done.size() << endl;*/
    }
    void AddCustomer(ACustomer c){
        // zamknu mutexem a pridam zakaznika do fronty nebo vectoru
        //cout << "AddCustomer running" << endl;
        //sem_wait(&sem_free);
        Customer* cust = new Customer(c);
        customers.push_back(cust);
        
        producents.emplace_back(&CRig::GivingThreadFIT, this, ref(*cust));
        producents.emplace_back(&CRig::GivingThreadCVUT, this, ref(*cust));
        consumers.emplace_back(&CRig::GettingThread, this, ref(*cust));
        //sem_post(&sem_free);
    }
  private:
    // predava problem tezeni pracovnim vlaknum
    void GivingThreadFIT(Customer& c){
        //cout << "GivingThreadFIT: Start" << endl;
        // vlozi problem do bufferu a nejspis pokud bude uspane working thread, tak ho probudi
        //int cnt = 0;
        AFITCoin x;
        while (x = c.GetCustomer().get()->FITCoinGen()) {
            /*if (cnt++ == 2)
                this_thread::sleep_for(chrono::seconds(1));*/
            sem_wait(&sem_free);
            {
                unique_lock<mutex> lock(mtx_buffer);
                //cout << "GivingThreadFIT: push_back" << endl;
                buffer.push_back(new Coin(ref(c), x));
            }
            {
                unique_lock<mutex> lock(c.mtx_cnt);
                c.request_cnt++;
            }
            //this_thread::sleep_for(chrono::seconds(1));
            sem_post(&sem_full);
        }
        /*sem_wait(&sem_free);
        {
            unique_lock<mutex> lock(mtx_buffer);
            cout << "GivingThreadFIT: push_back" << endl;
            buffer.push_back(new Coin(ref(c), x));
        }
        {
            unique_lock<mutex> lock2(c.mtx_cnt);
            c.request_cnt++;
        }*/
        {
            unique_lock<mutex> lock(c.mtx_fit_done);
            c.fit_done = true;
        }
        {
            unique_lock<mutex> lock(c.mtx_cnt);
            if (!c.request_cnt){
                lock.unlock();
                sem_post(&sem_full);
                sem_post(&c.sem_fill);
            }
        }
        // ceka na probuzeni
        //cout << "GivingThreadFIT: End" << endl;
    }
    // predava problem tezeni pracovnim vlaknum
    void GivingThreadCVUT(Customer& c){
        //cout << "GivingThreadCVUT: Start" << endl;
        // vlozi problem do bufferu a nejspis pokud bude uspane working thread, tak ho probudi
        //int cnt = 0;
        ACVUTCoin x;
        while (x = c.GetCustomer().get()->CVUTCoinGen()) {
            /*if (cnt++ == 2)
                this_thread::sleep_for(chrono::seconds(1));*/
            sem_wait(&sem_free);
            {
                unique_lock<mutex> lock(mtx_buffer);
                //cout << "GivingThreadCVUT: push_back" << endl;
                buffer.push_back(new Coin(ref(c), x));
            }
            {
                unique_lock<mutex> lock2(c.mtx_cnt);
                c.request_cnt++;
            }
            sem_post(&sem_full);
        }
        /*sem_wait(&sem_free);
        {
            unique_lock<mutex> lock(mtx_buffer);
            cout << "GivingThreadCVUT: push_back" << endl;
            buffer.push_back(new Coin(ref(c), x));
        }
        {
            unique_lock<mutex> lock(c.mtx_cnt);
            c.request_cnt++;
        }*/
        {
            unique_lock<mutex> lock(c.mtx_cvut_done);
            c.cvut_done = true;
        }        
        {
            unique_lock<mutex> lock(c.mtx_cnt);
            if (!c.request_cnt){
                lock.unlock();
                sem_post(&sem_full);
                sem_post(&c.sem_fill);
            }
        }
        // ceka na probuzeni
        //cout << "GivingThreadCVUT: End" << endl;
    }
    // ceka na vyreseni problemu a prebiraho zpet
    void GettingThread(Customer& c){
        //cout << "GettingThread: Start" << endl;
        Coin* coin;
        //int cnt = 0;
        while ( true ){
            /*if (cnt++ == 5)
                this_thread::sleep_for(chrono::seconds(1));*/
            sem_wait(&c.sem_fill);
            {
                unique_lock<mutex> lock(c.mtx_cnt);
                unique_lock<mutex> lock2(c.mtx_fit_done);
                unique_lock<mutex> lock3(c.mtx_cvut_done);
                //cout << "request_cnt: " << c.request_cnt << ", fit_done: " << c.fit_done << ", cvut_done: " << c.cvut_done << endl;
                if (!c.request_cnt && c.fit_done && c.cvut_done)
                    break;
            }
            {
                unique_lock<mutex> lock(c.mtx_cnt);
                unique_lock<mutex> lock2(c.mtx_fit_done);
                unique_lock<mutex> lock3(c.mtx_cvut_done);
                if(!c.request_cnt && (c.fit_done || c.cvut_done))
                    continue;
            }
            /*{
                unique_lock<mutex> lock(c.mtx_done);
                if (c.done.empty()){
                    lock.unlock();
                    sem_wait(&c.sem_fill);
                }
            }*/
            //cout << "GettingThread: Wake up" << endl;
            //this_thread::sleep_for(chrono::seconds(1));
            {
                unique_lock<mutex> lock(c.mtx_done);
                coin = c.done.front();
                c.done.erase(c.done.begin());
            }
            {
                unique_lock<mutex> lock2(c.mtx_cnt);
                c.request_cnt--;
            }
            sem_post(&c.sem_empty);
            //cout << "GettingThread: Accept" << endl;
            coin->Accept();
            delete coin;
            /*{
                unique_lock<mutex> lock(c.mtx_cnt);
                unique_lock<mutex> lock2(c.mtx_fit_done);
                unique_lock<mutex> lock3(c.mtx_cvut_done);
                if (!c.request_cnt && c.fit_done && c.cvut_done)
                    break;
            }*/
        }     
        //cout << "GettingThread: End" << endl;
    }
    void WorkingThread(void){
        //cout << "WorkingThread: Start" << endl;
        // bude si sahat pro problemy do fronty a volat nad nimi Solve(x)
        // !buffer.empty() || !stop_flag
        while (true) {
            sem_wait(&sem_full);
            {
                unique_lock<mutex> lock(mtx_buffer);
                unique_lock<mutex> lock2(mtx_stop);
                if (buffer.empty() && stop_flag)
                    break;
            }
            {
                unique_lock<mutex> lock(mtx_buffer);
                if (buffer.empty())
                    continue;
            }
            //sem_wait(&sem_full);
            //this_thread::sleep_for(chrono::seconds(1));
            //cout << "WorkingThread: In Cycle" << endl;
            Coin* c;
            {
                unique_lock<mutex> lock(mtx_buffer);
                c = buffer.front();
                buffer.erase(buffer.begin());
            }
            sem_post(&sem_free);
            //cout << "WorkingThread: Solve" << endl;
            c->Solve();

            sem_wait(&c->GetCustomer().sem_empty);
            {
                unique_lock<mutex> lock(c->GetCustomer().mtx_done);
                c->GetCustomer().done.push_back(c);
            }
            sem_post(&c->GetCustomer().sem_fill);
        }
        //cout << "WorkingThread: End" << endl;
    }
    // convert the vector<uint8_t> to string
    static void ConvertToString(const vector<uint8_t>& vec, string& str){
        // *it je typu uint8_t
        for (auto it = vec.end() - 1; it != vec.begin() - 1; --it){
            for (int i = 0; i < 8; ++i){
                str += ((*it & (1 << (7 - i))) ? "1" : "0");
            }
        }
    }
    // generate cartesian multiple of prefix and sufix
    static uint64_t GeneratePrefixSufix(const string& str, const int& min, const int& max){
        uint64_t cnt = 0;
        int res = 0;
        for (unsigned int i = 0; i < str.size(); ++i){
            for (unsigned int j = str.size(); j > 0; --j){
                res = lev_dist(str.substr(0, i + 1), str.substr(j - 1, str.size()));
                if (min <= res && res <= max)
                    cnt++;
            }
        }
        return cnt;
    }
    // method with levenshtein distance algorithm
    // source: https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C++
    static int lev_dist(const string& str1, const string& str2){
        const size_t len1 = str1.size(), len2 = str2.size();
        std::vector<std::vector<unsigned int>> dp(len1 + 1, std::vector<unsigned int>(len2 + 1));

        dp[0][0] = 0;
        for (unsigned int i = 1; i <= len1; i++) 
            dp[i][0] = i;
        for (unsigned int i = 1; i <= len2; i++)
            dp[0][i] = i;

        for (unsigned int i = 1; i <= len1; ++i){
            for (unsigned int j = 1; j <= len2; ++j){
                dp[i][j] = min({ dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + (str1[i - 1] == str2[j - 1] ? 0 : 1) });
            }
        }

        return dp[len1][len2];
    }
    // vraci pocet fixnich bitu
    static uint32_t SimplifiedVector(const vector<uint32_t>& vect, vector<uint32_t>& vec){
        array<bool, 32> fix_array; 
        fix_array.fill(true);
        int cnt = 0; 
        // kdyz je ve vektoru vice polozek, tak ho projdu a zjistim pocet fixnich bitu a zaznamenam si je do pole
        if (vect.size() > 1){
            bool first_step = true;
            uint32_t first = vect.front();
            uint32_t res = 0;
            for (const auto& v : vect){
                if(first_step){
                    first_step = false;
                    continue;
                }
                res = first ^ v;
                for (int i = 0; i < 32; i++){
                    if (res & (1 << i)){
                        fix_array[i] = false;
                    }
                }
            }            
            // zde provadim posun bitu smerem doprava, aby ve vysledku mohl v konecnem cyklu porovnavat mene hodnot
            for (const auto& v : vect){
                uint32_t tmp = 0;
                cnt = 0;
                for (int i = 0; i < 32; i++){
                    if (!fix_array[i] && (v & (1 << i))){
                        tmp += (0x1 << cnt);
                        cnt++;
                    } else if (!fix_array[i] && !(v & (1 << i))){
                        cnt++;
                    }
                }
                vec.push_back(tmp);
            }
        } else {
            fix_array.fill(false);
            uint32_t v = vect.front();
            uint32_t tmp = 0;
            for (int i = 0; i < 32; i++){
                if (!fix_array[i] && (v & (1 << i))){
                    tmp += (0x1 << cnt);
                    cnt++;
                }
            }
            vec.push_back(tmp);
        }
        return 32 - cnt;
    }
    static int FindDistance(uint32_t v, uint32_t i){
        int distance = 0;
        uint32_t res = i ^ v;
        for (uint32_t j = 0; j < 32; j++){
            if (res & (1 << j))
                distance++;
        }
        return distance;
    }
    static uint64_t CombinatoricNumber(uint32_t n, uint32_t k){
        if (n < k)
            return 0;        
        if (n < 2 * k)
            k = n - k;
        if (k == 0)
            return 1;

        uint64_t res = n;
        for (uint32_t i = 2; i <= k; ++i){
            res = (res * (n - i + 1)) / i;
        }
        return res;
    }
    static uint64_t FindDist(const vector<uint32_t>& vec, const int& distMax, const uint32_t& fix){
        uint64_t cnt = 0;
        for (uint32_t i = 0; i < pow(2, 32 - fix); i++){
            int dist = 0;
            for (const auto& v : vec){
                dist = max(dist, FindDistance(v, i));
            }
            if (dist == distMax)
                cnt++;
            else if (dist < distMax){
                for (int j = 0; j <= (distMax - dist); j++){   // fix = pocet fixnich bitu
                    cnt += CombinatoricNumber(32 - (32 - fix), j);
                }
            }
        }
        return cnt;
    }

    class Coin;
    class Customer;
    
    class Customer{
    public:
        Customer(void): request_cnt(0), fit_done(false), cvut_done(false) {
            sem_init(&sem_empty, 0, BUFFER_SIZE_CUSTOMER);
            sem_init(&sem_fill, 0, 0);
        }
        Customer(const ACustomer& x): Customer() {
            c = x;
        }
        Customer(const Customer& x) {
            c = x.c;
            request_cnt = x.request_cnt;
            fit_done = x.fit_done;
            cvut_done = x.cvut_done;
            done = x.done;
            sem_fill = x.sem_fill;
            /*sem_destroy(&sem_fill);
            sem_init(&sem_fill, 0, 0);*/
        }
        ~Customer() {
            sem_destroy(&sem_empty);
            sem_destroy(&sem_fill);
        }
        Customer& operator=(const Customer& o){
            c = o.c;
            request_cnt = o.request_cnt;
            fit_done = o.fit_done;
            cvut_done = o.cvut_done;
            done = o.done;
            sem_fill = o.sem_fill;
            return *this;
        }
        ACustomer GetCustomer() {
            return c;
        }
        vector<Coin*>    done;
        int              request_cnt;
        bool             fit_done, cvut_done;
        mutex            mtx_done, mtx_cnt, mtx_fit_done, mtx_cvut_done;
        sem_t            sem_fill, sem_empty;
    private:
        ACustomer       c;
    };
    class Coin{
    public:
        Coin(void) { }
        Coin(Customer& c, AFITCoin x) { 
            fit = true;
            f = x;
            this->customer = &c;
        }
        Coin(Customer& c, ACVUTCoin x) {
            fit = false;
            this->c = x;
            this->customer = &c;
        }
        void Solve(void){
            //cout << "fit: " << fit << ", fit_done: " << fit_done << ", cvut_done: " << cvut_done << endl;
            if (fit){
                if (!f){
                    return;
                }
                CRig::Solve(f);
            } else {
                if (!this->c){
                    return;
                }
                CRig::Solve(c);
            }
        }
        void Accept(void){
            if (fit){
                if (!f)
                    return;
                customer->GetCustomer()->FITCoinAccept(f);
            } else {
                if (!this->c)
                    return;
                customer->GetCustomer()->CVUTCoinAccept(c);
            }
        }
        Coin& operator=(const Coin& o){
            this->customer = o.customer;
            this->f = o.f;
            this->c = o.c;
            this->fit = o.fit;

            return *this;
        }
        Customer& GetCustomer(void){
            return *customer;
        }
    private:
        Customer*    customer;
        AFITCoin     f;
        ACVUTCoin    c;
        bool         fit;
    };
    
    int thrCnt;
    // buffer added customers waiting for perform
    vector<Customer*>    customers;
    // buffers
    vector<Coin*>        buffer;
    // vector of all threads
    vector<thread>      producents;
    vector<thread>      consumers;
    vector<thread>      workers;
    mutex               mtx_buffer, mtx_stop;
    sem_t               sem_free, sem_full;
    bool                stop_flag;
};


#ifndef __PROGTEST__
#include "test.inc"
#endif /* __PROGTEST__ */
