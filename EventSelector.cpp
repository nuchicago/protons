//EventSelector class



#include <iostream>



class EventSelector{

	

	public:
		EventSelector ();
		EventSelector (int);
		int EventNtuple;
		int PrintRepr ();
		



};


EventSelector::EventSelector(int n){
	EventNtuple = n;
	std::cout << "Ntuple assigned \n";
	

}

EventSelector::EventSelector(){
    std::cout << "No Ntuple given \n";
	
}


int EventSelector::PrintRepr(){

    std::cout << "EventSelector class \n";
	return 0;
}

int main(){

	EventSelector t1;
	t1.PrintRepr();
	return 0;



}