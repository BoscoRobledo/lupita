#ifndef SC_STOP_CRITERION_H_INCLUDED
#define SC_STOP_CRITERION_H_INCLUDED
class StopCriterion{
    public:
        virtual bool stopNow(int k) = false;
        int getMaxk()
        {
            return maxk;
        }
        void setMaxk(int k)
        {
            maxk=k;
        }
        virtual ~StopCriterion(){};
    	
    private:
        int maxk;
	
};
#endif
