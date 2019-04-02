#ifndef _LATTICE_H_
#define _LATTICE_H_
#include <vector>
#include <array>
#include <stdexcept>

        
class DiscreteLattice{
    public:
        //Arguments:
        //  dim: dimension of the lattice (defult 3)
        DiscreteLattice(unsigned int dim=3);
        ~DiscreteLattice();

        //Set position on the lattice
        //Arguments: 
        //  pos:    a std::vector retaining the position to be set. 
        //          Must be of the same dimension of the lattice, if wrong a
        //          std::range_error exception will be thrown.
        void setPosition(const std::vector<int>& pos);

        //Returns a std::vector containing the actual position on lattice
        std::vector<int> getPosition() const;
        
        //Move the actual position among the ith axis (with numeration starting
        //from 1, e.g: for a 3D lattice axis would be 1,2,3), forward if second
        //arg is true, backward elsewhere. 
        //
        //Arguments:
        //  axis:       ith axis to move on, starting from 1
        //  forward:    move forward if true, backward if false
        void move(unsigned int axis, bool forward);

        //Return squared distance of position from the origin
        int rsquared() const; 

        //move position to origin
        void origin();
    private:
        std::vector<int> position;
};

class Dense3DRW{
    public:
        //Dense "lattice" 3-dimensional random walks.
        //
        //step_size is the a lenght of the lattice. Each time the walker moves 
        //in the given direction with a step of length step_size
        //Arguments:
        //  step_size:  length of moves
        Dense3DRW(double step_size=1);
        ~Dense3DRW();
         //Set position on the lattice
        //Arguments: 
        //  pos:    a std::array retaining the position to be set. 
        void setPosition(const std::array<double,3>& pos);

        //Returns a std::array containing the actual position on lattice
        std::array<double,3> getPosition() const;
        
        //Moves the walker in the direction given by polar angles theta and phi
        //by a distance of step_size (given in the constructor of this class).
        //Arguments:
        //  theta:      tetha angle, to be in range [0,pi]
        //  phi:        phi angle, to be in range [0, 2pi]
        void move(double theta, double phi);

        //Return squared distance of position from the origin
        double rsquared() const; 

        //move position to origin
        void origin();
    private:
        std::array<double,3> position;
        double stepsize;
};

        



inline std::vector<int> DiscreteLattice::getPosition() const{
    return position;
}

inline std::array<double,3> Dense3DRW::getPosition() const{
    return position;
}
#endif //_LATTICE_H_
