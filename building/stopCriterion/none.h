#ifndef SC_NONE_H_INCLUDED
#define SC_NONE_H_INCLUDED
class None: public StopCriterion{
    public:
        bool stopNow(int k)
        {
            return false || k>maxk;
        }
        virtual ~None(){};	
};
#endif
