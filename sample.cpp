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

class CRig
{
  public:
    static void Solve(ACVUTCoin x){
        string str = "";
        vector<uint8_t>& vec = x.get()->m_Data;
        ConvertToString(vec, str);
        // pouze pro debugovani
        //cout << str << endl;
        x.get()->m_Count = GeneratePrefixSufix(str, x.get()->m_DistMin, x.get()->m_DistMax);
    }
    static void Solve(AFITCoin x){
        if (x.get()->m_Vectors.size() == 0){
            return;
        }
        vector<uint32_t> vec;
        uint32_t fix = 0;
        // SimplifiedVector
        fix = SimplifiedVector(x.get()->m_Vectors, vec);
        x.get()->m_Count = FindDist(vec, x.get()->m_DistMax, fix);
    }
    CRig(void){

    }
    ~CRig(void){

    }
    void Start(int thrCnt){

    }
    void Stop(void){
        // joinovat vlakna nebo je ?? killnout ??
    }
    void AddCustomer(ACustomer c){
        // zamknu mutexem a pridam zakaznika do fronty nebo vectoru
    }
  private:
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
    // todo
    // add vector for custumers
    // add mutex
};


#ifndef __PROGTEST__
#include "test.inc"
#endif /* __PROGTEST__ */
