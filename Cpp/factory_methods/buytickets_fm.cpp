//modified from https://www.geeksforgeeks.org/design-patterns-set-2-factory-method/

// C++ program to demonstrate factory method design pattern
#include <iostream>
using namespace std;

enum TicketType {
	FlightTicket, BusTicket, TrainTicket
};

// Library classes
class Tickets {
public:
	virtual void printItinerary() = 0;
  virtual void name(std::string&) = 0;
  virtual void origin(std::string&) = 0;
  virtual void destination(std::string&) = 0
	static Tickets* Create(TicketType type);
};

class FlightTicket : public Tickets {
public:
	void printVehicle() {
		cout << "I am two wheeler" << endl;
	}
};
class BusTicket : public Tickets {
public:
	void printVehicle() {
		cout << "I am three wheeler" << endl;
	}
};
class TrainTicket : public Tickets {
	public:
	void printVehicle() {
		cout << "I am four wheeler" << endl;
	}
};

// Factory method to create objects of different types.
// Change is required only in this function to create a new object type
Vehicle* Vehicle::Create(VehicleType type) {
	if (type == VT_TwoWheeler)
		return new TwoWheeler();
	else if (type == VT_ThreeWheeler)
		return new ThreeWheeler();
	else if (type == VT_FourWheeler)
		return new FourWheeler();
	else return NULL;
}

// Client class
class Client {
public:

	// Client doesn't explicitly create objects
	// but passes type to factory method "Create()"
	Client()
	{
		VehicleType type = VT_FourWheeler;
		pVehicle = Vehicle::Create(type);
	}
	~Client() {
		if (pVehicle) {
			delete[] pVehicle;
			pVehicle = NULL;
		}
	}
	Vehicle* getVehicle() {
		return pVehicle;
	}

private:
	Vehicle *pVehicle;
};

// Driver program
int main() {
	Client *pClient = new Client();
	Vehicle * pVehicle = pClient->getVehicle();
	pVehicle->printVehicle();
	return 0;
}


