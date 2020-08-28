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

// to continue to decrypt the intermediate monomials
int SecondBackExpandPolynomial( int rounds, bitset<288> final, vector<bitset<288> > & term )
{
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2);//focus on finding best solutions 
    env.set(GRB_IntParam_PoolSolutions, MAX); // try to find 2000000
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

    for ( int i = 0; i < 288; i++ )
        if ( final[i] == 0 )
            model.addConstr( works[i] == 0 );
        else
            model.addConstr( works[i] == 1 );

    model.update();
    model.optimize();

    map<bitset<288>, int, cmp288> counterMap; 

    if( model.get( GRB_IntAttr_Status ) == GRB_OPTIMAL )
    {
        double time = model.get(GRB_DoubleAttr_Runtime );
        //cout << "Time Used: " << time << "sec" << endl;
        
        int solCount = model.get(GRB_IntAttr_SolCount);
        //cout << "Raw Solutions: " << solCount << endl;

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
    //cout << "Exact terms: " << term.size() << endl;
}

int  MidSolutionCounter( int rounds, const bitset<285> & start, const bitset<288> & last, double& time )
{
    //setting the enviroment
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, 48);
    env.set(GRB_IntParam_PoolSearchMode, 2);//finding n best solutions 
    env.set(GRB_IntParam_PoolSolutions, MAX); // try to find 2000000
    GRBModel model = GRBModel(env);

    vector<GRBVar> s(288);
    for (int i = 0; i < 288; i++)
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);

    // set the initial state, note the last three bits are constant 1
    for ( int i = 0; i < 285; i++ )
        if ( start[i] == 0 )
            model.addConstr( s[i] == 0 );
        else
            model.addConstr( s[i] == 1 );

    //propagate the monomial
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

    // the tail is the monomial in U
    for ( int i = 0; i < 288; i++ )
        if ( last[i] == 1)
            model.addConstr( works[i] == 1 );
        else
            model.addConstr( works[i] == 0 );

    model.set(GRB_DoubleParam_TimeLimit, 300.0 );
    model.optimize();
    int solnum = 0;
    if ( model.get( GRB_IntAttr_Status ) == GRB_TIME_LIMIT )
    {
	    //cout << "-------------------------------------------------------------- EXPAND" << endl;
	    vector<bitset<288>> T;
	    int re = 0;
        do 
	    {
		    T.clear();
            re ++;
            SecondBackExpandPolynomial(re, last, T );
	    }while( T.size() <= 16 && ( re + 10 ) < rounds );
	    int tsize = T.size();
	    int c = 0;
        for ( auto it : T )
        {
            //cout << c << " out of " << tsize << "| Depth " << depth << endl;
            c++;
	        double mytime; 
            solnum += MidSolutionCounter( rounds - re, start, it, mytime );    
        }
	    depth --;
    }
    else 
    {
         time = model.get(GRB_DoubleAttr_Runtime );
         int solCount = model.get(GRB_IntAttr_SolCount);

         return solCount;
    }
}

/*
 * rounds: the round we want to decrypt
 * term: all terms in U
 */
void BackExpandPolynomial( int rounds, vector<bitset<288> > & term )
{
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2); // PoolSearchMode = 2 is better
    env.set(GRB_IntParam_PoolSolutions, MAX); // set PoolSolution Size
    GRBModel model = GRBModel(env);

    vector<GRBVar> s(288);
    for (int i = 0; i < 288; i++)
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);

    // propagate the monomials
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

    // output function
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
        //cout << "Time Used: " << time << "sec" << endl;
        
        int solCount = model.get(GRB_IntAttr_SolCount);
        //cout << "Raw Solutions: " << solCount << endl;

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
    //cout << "Exact terms: " << term.size() << endl;
}

/*
 * rounds: the target round
 * backrounds: decryption rounds 
 * Term: the monomials of trivium after decrypting backrounds Trivium
 * start: as a return value, return the confirmed monomial
 */
int UpBound( int rounds, int backrounds, vector<bitset<288>> & Term, bitset<285> & start )
{
    // for C++ API, the model is always declared as this, see Gurobi Manual 
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0); // close output to screen
    env.set(GRB_IntParam_Threads, 48 ); // threads = 48
    GRBModel model = GRBModel(env);

    vector<GRBVar> s(288); // the initial register
    for (int i = 0; i < 288; i++)
        s[i] = model.addVar(0, 1, 0, GRB_BINARY); // declare binary variables

    for ( int i = 0; i < 80; i++ ) // keys are all secret variables
        ;
    for ( int i = 80; i < 93; i++ )
        model.addConstr( s[i] == 0 ); // constant 0

    for ( int i = 93; i < 93 + 80; i++ ) // IVs are all public variables 
        ;

    for ( int i = 80 + 93; i < 285; i++ ) // constant 0
        model.addConstr( s[i] == 0 );

    for ( int i = 285; i < 285; i++ ) // constant 1
        ;

    // propagates the monomial, similar to HLM+20 
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

    // output function of trivium
    GRBLinExpr nk = 0;
    for ( int i = 0; i < 288; i++ )
        if ( (i == 65) || (i == 92) || (i == 161) || (i == 176) || (i == 242) || (i == 287))
            nk += works[i];
        else 
            model.addConstr( works[i] == 0);
    model.addConstr( nk == 1 );

    // the degree 
    GRBLinExpr nv = 0;
    for ( int i = 93; i < 93 + 80; i++ )
        nv += s[i];

    // max the degree as the objective function
    model.setObjective(  nv, GRB_MAXIMIZE );

    int size = Term.size();

    while ( true )
    {
        model.update();
        model.optimize();

        if ( model.get( GRB_IntAttr_Status ) == GRB_OPTIMAL )
        {
            double time = model.get(GRB_DoubleAttr_Runtime );
            // find a monomial x^u \rightsquigarrow z 
            //cout << "Candidate Found: " 
            //     <<  model.getObjective().getValue() 
            //     <<  " |  Time Used: " 
            //     << time << "sec" << endl;
        
            for ( int j = 0; j < 285; j++ ) 
                if ( round( s[j].get( GRB_DoubleAttr_Xn ) ) == 1 )  
                    start[j] = 1;
                else 
                    start[j] = 0;

            int solCount = 0;
            //cout << "We need to count " << size <<  " terms" << endl;

            // divide-and-conquer
            for ( int i = 0; i < size; i++ )
            {
                int solnum = MidSolutionCounter( rounds - backrounds, start, Term[i], time );
                solCount += solnum;
                //cout << "Term: " << i << " Solnum: " << solnum << " Total: " << solCount  << " |  Time: "  << time << endl;
            }

            //confirmed
            if ( solCount % 2 == 1 )
                return round ( model.getObjective().getValue() );
            else  // eliminating the solution
            {
                GRBLinExpr excludeCon = 0;
                for (int i = 0; i < 285; i++)
                    if (start[i] == 1)
                        excludeCon += (1 - s[i]);
                    else
                        excludeCon += s[i];
                model.addConstr( excludeCon >= 1);
            }
        }
    }
}

/*
 * Store the degree and the confirmed monomial
 */
struct DV
{
    int deg;
    bitset<285> vec;
};

int main()
{
    DV DEG[900];
    int MID = 1; // For divide-and-conquer
    int str = 2, end = 834; // from round 2 - 834

    vector<bitset<288>> Term;
    BackExpandPolynomial( MID, Term ); // Expand trivium
    
    for ( int i = str; i < end; i++ )
    {
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Round: " << i << endl;
        cout << getCurrentSystemTime() << endl;
        DEG[i].deg = UpBound( i, MID, Term, DEG[i].vec );
        cout << "-----------------------------------------------------------------------" << endl;
        cout << "Round: " << i << " | "  << "Upbound: " <<  DEG[i].deg << endl;; 
	    //cout << DEG[i].vec << endl;
        for ( int j = 0; j < 80; j++ )
            if ( DEG[i].vec[j] == 1 )
                cout << "k" << j << " ";

        cout << "  |  ";

        for ( int j = 0; j < 80; j++ )
            if ( DEG[i].vec[j + 93] == 1 )
                cout << "x" << j << " ";
        cout << endl;

        cout << "-----------------------------------------------------------------------" << endl;
        cout << getCurrentSystemTime() << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    }

    ofstream deglog;
    deglog.open( "deg.txt" );
    for ( int i = str; i < end; i++ )
    {
        deglog << "Round: "  << i << " |   Degree: " << DEG[i].deg << endl;
        for ( int j = 0; j < 80; j++ )
            if ( DEG[i].vec[j] == 1 )
                deglog << "k" << j << " ";

        deglog << "  |  ";

        for ( int j = 0; j < 80; j++ )
            if ( DEG[i].vec[j + 93] == 1 )
                deglog << "x" << j << " ";
        deglog << endl;
    }
}

