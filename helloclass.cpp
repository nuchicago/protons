//Hello Class!



#include <iostream>


class sample_class{



	public:
		int PrintRepr ();

		




};
int sample_class::PrintRepr(){

    std::cout << "Hello Class Template \n";
	return 0;
}

int main(){

	sample_class t1;
	t1.PrintRepr();
	return 0;



}