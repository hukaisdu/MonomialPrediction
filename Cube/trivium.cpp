#include<iostream>
#include<cstdio>
#include<bitset>
#include<vector>
#include<set>
#include<map>
#include<cmath>
#include<fstream>
#include<chrono>
#include"gurobi_c++.h" 

using namespace std;

#define THREAD 28

int depth = 0;

const int MAX = 200000000; // the maximum value of PoolSearchMode, P625

string getCurrentSystemTime()
{
    auto tt = chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    struct tm* ptm = localtime(&tt);
    char date[60] = { 0 };
    sprintf(date, "%d-%02d-%02d-%02d:%02d:%02d", (int)ptm->tm_year + 1900, (int)ptm->tm_mon + 1, (int)ptm->tm_mday,
                                        (int)ptm->tm_hour, (int)ptm->tm_min, (int)ptm->tm_sec);
    return string(date);
}

struct cmp288
{
    bool operator()( const bitset<288> & a, const bitset<288> & b ) const
    {
        for ( int i = 0; i < 288; i++ )
            if ( a[i] < b[i] ) return true;
            else if ( a[i] > b[i] ) return false;
        return false; // equal
    }
};

struct cmp285
{
    bool operator()( const bitset<285> & a, const bitset<285> & b ) const
    {
        for ( int i = 0; i < 285; i++ )
            if ( a[i] < b[i] ) return true;
            else if ( a[i] > b[i] ) return false;
        return false; // equal
    }
};

void triviumCore(GRBModel& model, vector<GRBVar>& x, int i1, int i5, int i2, int i3, int i4)
{    
   int Ineq[][11] = {
 {0, -1, -1, 0, -1, -1, 1, 1, 0, 1, 1},
 {0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0},
 {0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0},
 {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1},
 {0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0},
 {0, -1, 0, -1, -1, -1, 1, 0, 1, 1, 1},
 {0, 0, -1, 1, 0, 0, 0, 1, 0, 0, 0},
 {0, 0, 1, -1, 0, 0, 0, 0, 1, 0, 0},
 {2, 0, 1, 0, 1, 0, -1, 0, 0, -1, -1},
 {0, 1, 1, 0, 1, 1, 0, 0, 0, -1, 0},
 {3, 1, 0, 0, 1, 0, 0, -1, -1, -1, -1},
 {2, 0, 0, 1, 1, 0, -1, 0, 0, -1, -1},
 {0, 1, 0, 1, 1, 1, 0, 0, 0, -1, 0},
 {3, 0, 0, 0, 1, 1, -1, -1, -1, -1, 0}
 };  
 GRBVar y1 = model.addVar(0, 1, 0, GRB_BINARY);
 GRBVar y2 = model.addVar(0, 1, 0, GRB_BINARY);
 GRBVar y3 = model.addVar(0, 1, 0, GRB_BINARY);
 GRBVar y4 = model.addVar(0, 1, 0, GRB_BINARY);
 GRBVar y5 = model.addVar(0, 1, 0, GRB_BINARY);
    
 for ( auto it : Ineq )
    model.addConstr( it[0] + it[1] * x[i1] + it[2] * x[i2] + it[3] * x[i3] + it[4] * x[i4] + it[5] * x[i5] + it[6] * y1 + it[7] * y2 + it[8] * y3 + it[9] * y4 +  it[10] * y5 >= 0 );

 x[i1] = y1;
 x[i2] = y2;
 x[i3] = y3;
 x[i4] = y4;
 x[i5] = y5;
}

/*
void triviumCore(GRBModel& model, vector<GRBVar>& x, int i1, int i2, int i3, int i4, int i5)
{
    GRBVar y1 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar y2 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar y3 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar y4 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar y5 = model.addVar(0, 1, 0, GRB_BINARY);

    GRBVar z1 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar z2 = model.addVar(0, 1, 0, GRB_BINARY);

    // z3 and z4 are not needed, since z3 = z4 = a
    GRBVar a = model.addVar(0, 1, 0, GRB_BINARY);

    //copy
    model.addConstr(y1 <= x[i1]);
    model.addConstr(z1 <= x[i1]);
    model.addConstr(y1 + z1 >= x[i1]);

    //copy
    model.addConstr(y2 <= x[i2]);
    model.addConstr(z2 <= x[i2]);
    model.addConstr(y2 + z2 >= x[i2]);

    //copy
    model.addConstr(y3 <= x[i3]);
    model.addConstr(a <= x[i3]);
    model.addConstr(y3 + a >= x[i3]);
    
    //copy
    model.addConstr(y4 <= x[i4]);
    model.addConstr(a <= x[i4]);
    model.addConstr(y4 + a >= x[i4]);
    //XOR
    model.addConstr(y5 == x[i5] + a + z1 + z2);

    x[i1] = y1;
    x[i2] = y2;
    x[i3] = y3;
    x[i4] = y4;
    x[i5] = y5;
}
*/

int  MidSolutionCounter( int rounds, bitset<80> cube, const bitset<288> & last, 
map<bitset<285>, int, cmp285> & counterMap, ostream & f = cout )
{
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, THREAD);
    env.set(GRB_IntParam_PoolSearchMode, 2);//focus on finding n best solutions 
    env.set(GRB_IntParam_MIPFocus, 3 );
    env.set(GRB_IntParam_PoolSolutions, MAX); // try to find 2000000
    GRBModel model = GRBModel(env);

    vector<GRBVar> s(288);
    for (int i = 0; i < 288; i++)
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);

    for ( int i = 0; i < 80; i++ )
        if ( cube[i] == 0 )
            model.addConstr( s[i + 93] == 0 );
        else
            model.addConstr( s[i + 93] == 1 );

    // key, last three bits
    for ( int i = 80; i < 93; i++ )
        model.addConstr( s[i] == 0 );
    for ( int i = 93 + 80; i < 285; i++ )
        model.addConstr( s[i] == 0 );

    vector<GRBVar> works = s;
    for (int r = 0; r < rounds; r++) 
    {
        triviumCore(model, works, 65, 170, 90, 91, 92);
        triviumCore(model, works, 161, 263, 174, 175, 176);
        triviumCore(model, works, 242, 68, 285, 286, 287);
            
        vector<GRBVar> temp = works;
        for (int i = 0; i < 288; i++) 
            works[(i + 1) % 288] = temp[i];
    }

    // Output function
    for ( int i = 0; i < 288; i++ )
        if ( last[i] == 1)
            model.addConstr( works[i] == 1 );
        else
            model.addConstr( works[i] == 0 );

    GRBLinExpr nk = 0;
    for ( int i = 0; i < 80; i++ )
        nk += s[i];
    model.setObjective( nk, GRB_MAXIMIZE );

    cout << getCurrentSystemTime() << endl;
    //f << getCurrentSystemTime() << endl;
    model.optimize();

    double time = model.get(GRB_DoubleAttr_Runtime );
    cout << "Rounds: " << rounds << "  Time Used: " << time << "sec" << endl;
    //f << "Rounds: " << rounds << "  Time Used: " << time << "sec" << endl;
    
    int solCount = model.get(GRB_IntAttr_SolCount);
    cout << "---------------------------------------Raw Solutions: " << solCount << endl;
    //f << "----------------------------------------Raw Solutions: " << solCount << endl;

    bitset<285> start;
    for ( int i = 0; i < solCount; i++ )
    {
        model.set(GRB_IntParam_SolutionNumber, i );
        for ( int j = 0; j < 285; j++ ) 
            if ( round( s[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                start[j] = 1;
            else 
                start[j] = 0;
        counterMap[start]++;
    }
    return 0;
}

int BackExpandPolynomial( int rounds, vector<bitset<288> > & term )
{
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2); 
    env.set(GRB_IntParam_PoolSolutions, MAX); 

    GRBModel model = GRBModel(env);

    vector<GRBVar> s(288);
    for (int i = 0; i < 288; i++)
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);

    vector<GRBVar> works = s;
    for (int r = 0; r < rounds; r++) 
    {
        triviumCore(model, works, 65, 170, 90, 91, 92);
        triviumCore(model, works, 161, 263, 174, 175, 176);
        triviumCore(model, works, 242, 68, 285, 286, 287);
            
        vector<GRBVar> temp = works;
        for (int i = 0; i < 288; i++) 
            works[(i + 1) % 288] = temp[i];
    }

    GRBLinExpr nk = 0;
    for ( int i = 0; i < 288; i++ )
        if ( (i == 65) || (i == 92) || (i == 161) || (i == 176) || (i == 242) || (i == 287))
            nk += works[i];
        else 
            model.addConstr( works[i] == 0);
    model.addConstr( nk == 1 );

    model.update();
    model.optimize();

    map<bitset<288>, int, cmp288> counterMap; 

    if( model.get( GRB_IntAttr_Status ) == GRB_OPTIMAL )
    {
        double time = model.get(GRB_DoubleAttr_Runtime );
        cout << "Time Used: " << time << "sec" << endl;
        
        int solCount = model.get(GRB_IntAttr_SolCount);
        cout << "Raw Solutions: " << solCount << endl;

        bitset<288> start;
        for ( int i = 0; i < solCount; i++ )
        {
            model.set(GRB_IntParam_SolutionNumber, i );

            for ( int j = 0; j < 288; j++ ) 
                if ( round( s[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                    start[j] = 1;
                else 
                    start[j] = 0;
            counterMap[start]++;
        }
    }
    else if( model.get( GRB_IntAttr_Status ) == GRB_INFEASIBLE )
    {
        cout << "No terms " << endl;
        exit(-2);
    }
    else
    {
        cout << "Other status " << GRB_IntAttr_Status <<  endl;
        exit(-1);
    }

    for ( auto it : counterMap )
        if ( it.second % 2 == 1 )
            term.push_back( it.first );
    cout << "Exact Solutions: "  << term.size() << endl;
    return 0;
}


int main( int argc, char * argv[] )
{
    if ( argc != 3 )
    {
        cout << "Usage: ./trivium round cube_index" << endl;
        cout << endl;
        cout << "Valid Parameter: " << endl
             << "round = 840" << endl
             << "cube_index = 1: [0,1,...,79]/{70, 72, 74, 76, 78}" << endl
             << "cube_index = 2: [0,1,...,79]/{72, 74, 76, 78}" << endl
             << "cube_index = 3: [0,1,...,79]/{70, 74, 76, 78}" << endl
             << endl
             << "round = 841" << endl
             << "Cube_index = 1: [0,1,...,79]/{70, 72, 76, 78}" << endl
             << "Cube_index = 2: [0,1,...,79]/{72, 76, 78}" << endl
             << endl
             << "round = 842" << endl
             << "Cube_index = 1: [0,1,...,79]/{72, 74, 76, 78}" << endl
             << "Cube_index = 2: [0,1,...,79]/{74, 76, 78}" << endl;
        return -1;
    }
    int round = atoi( argv[1] );
    int index = atoi( argv[2] );
    int MID = 200;

    if ( !( round == 840 || round == 841 || round == 842 ) )
    {
        cout << "Usage: ./trivium round cube_index" << endl;
        cout << "The first parameter can only be 840, 841, 842 rounds" << endl;
        return -2;
    }
    else
    {
        if ( round == 840 )
        {
            if ( ! ( ( index == 1 ) || ( index == 2) || ( index == 3 ) ) ) 
            {
                cout << "Usage: ./trivium round cube_index" << endl;
                cout << "The second parameter can only be 1, 2, 3 for round = 840" << endl;
                return -3;
            }
        }
        if ( round == 841 )
        {
            if ( ! ( ( index == 1 ) || ( index == 2) ) ) 
            {
                cout << "Usage: ./trivium round cube_index" << endl;
                cout << "The second parameter can only be 1, 2 for round = 841" << endl;
                return -3;
            }
        }
        if ( round == 842 )
        {
            if ( ! ( ( index == 1 ) || ( index == 2) ) ) 
            {
                cout << "Usage: ./trivium round cube_index" << endl;
                cout << "The second parameter can only be 1, 2 for round = 842" << endl;
                return -3;
            }
        }
    }
    vector<int> I;
    if ( round == 840 )
    {
        MID = 200;
        if ( index == 1 )
        {
            int c[] = {70, 72, 74, 76, 78};
            for ( auto it : c )
                I.push_back( it );
        }
        if ( index == 2 )
        {
            int c[] = {72, 74, 76, 78};
            for ( auto it : c )
                I.push_back( it );

        }
        if ( index == 3 )
        {
            int c[] = {70, 74, 76, 78};
            for ( auto it : c )
                I.push_back( it );
        }
    }
    else  if ( round == 841 )
    {
        MID = 250;
        if ( index == 1 )
        {
            int c[] = {70, 72, 76, 78};
            for ( auto it : c )
                I.push_back( it );
        }
        if ( index == 2 )
        {
            static int c[] = {72, 76, 78};
            for ( auto it : c )
                I.push_back( it );
        }
    }
    else 
    {
        MID = 300;
        if ( index == 1 )
        {
            int c[] = {72, 74, 76, 78};
            for ( auto it : c )
                I.push_back( it );
        }
        if ( index == 2 )
        {
            static int c[] = {74, 76, 78};
            for ( auto it : c )
                I.push_back( it );
        }
    }

    int ROUND = round;

    // the start time
    auto start = chrono::steady_clock::now();

    vector<bitset<288>> midTerms;
    vector<bitset<80>> initialTerm;
    map< bitset<285>, int, cmp285 > gross, temp;

    bitset<80> cube;
    for ( int i = 0; i < 80; i++ )
        cube[i] = 1;
    for ( auto it : I )
        cube.reset(it); 

    cout << "Round: " << round << endl; 
    cout << "Cube: " << cube << endl;
    cout << "Thread: " << cube << endl;

    BackExpandPolynomial( MID, midTerms );

    int count = 0;
    for ( auto it : midTerms )
    {	
	    depth = 0;
        cout << "Trivium-" << ROUND << "   Term: " << count << "-th" << endl;
	    count++;
        temp.clear();
        MidSolutionCounter(ROUND - MID, cube, it, temp); 
        for ( auto jt : temp )
            gross[jt.first] += jt.second;
    }

    cout << "Raw Terms " << gross.size() << endl;

    for ( auto it: gross )
    {
        cout << it.first << " [" << it.second << "]" << endl;
    }

    for ( auto it : gross )
        if ( it.second % 2 == 1 )
        {
            bitset<80> tmp(0);
            for ( int i = 0; i < 80; i++ )
                tmp[i] = it.first[i];
            
            initialTerm.push_back( tmp );
        }

    
    ofstream res;
	res.open( "result.txt" );

    res  << "ROUND : " << ROUND << endl 
         << "Cube: " << cube << endl
         << "Terms: " << initialTerm.size() << endl;

    cout  << "ROUND : " << ROUND << endl 
         << "Cube: " << cube << endl
         << "Terms: " << initialTerm.size() << endl;

    for ( auto it : initialTerm )
    {
        cout << it << endl;
	    res << it << endl;
        int flg = 0;
        for ( int i = 0; i < 80; i++ )
        {
            if ( it[i] == 1 )
            {
                cout << "k" << i << ' ';
                flg = 1; // not all zero
            }
        }
        if ( flg == 0 ) // all are zero
            cout << "1"; 
        cout << endl;
    }

    auto end = chrono::steady_clock::now();
    auto time = chrono::duration<double> ( end - start );
    cout << getCurrentSystemTime() << endl;

    res << "Time: " << time.count() << " seconds" << endl;
    cout << "Time: " << time.count() << " seconds" << endl;
    res.close();
}

