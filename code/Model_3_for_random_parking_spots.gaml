/**
* Name: Distribution of advertising material
* Based on the internal empty template. 
* Author: Kilian Heyer
* 
* Purpose:
* This model was created during and for the Masther Thesis of Kilian Heyer.
* This version is for sampling random parking spots as first parking spot.
* Its purpose is to predict the time necessary to distibute advertising material (example: brochures) using nearest neighbor algorithm.
* This prediction does not correlate with reality very well and serves for comparison of different deliver areas with each other.
* 
* To run this model sucessfully you must give it the adresses, streets, start_point,road for deliverer as shp and the delivery_time as csv.
* An error can also appear if the files to save the data do not exist.    
* 
* 
*/


model Simulation_of_delivery_with_car


global {
    // loading the necessary files
    //the street_shapefile should contain connected lines, a column for the maximum speed in each direction with -1 in it if it is a one way street. It is used to create a graph for the streets for the vehicle. 
    //Example path: "D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502601/Straßennetzwerk_Lambert_31287_Gebiet_502601_für das vehicle.shp"
    file street_shapefile <- file("Path/to/your/shapefile.shp");
    
    //the adress_shapefile should contain points and aclomun for the number of households at the adress.
    //Example path: "D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502601/Adressen_502601.shp"
	file adress_shapefile <- file("Path/to/your/shapefile.shp");
	
	// the start_point should contain single point. 
	//Example path: "D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502601/FILIALEN_Salzburg.shp"
	file start_point <-file ("Path/to/your/shapefile.shp");
	
	//the building_shapefile is only cosmetic and not necessary. if included it should contain polygons.
	//Example path: "D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502601/Gebäude_502601.shp"
	file building_shapefile <-file ("Path/to/your/shapefile.shp");
	
	//this file should be a cas containing a header and in the first column continuous numbers from from 1 to the maximum number of households. In the second column should be the needed time for delivery ath the adress in minutes
	//Example path: "D:/Master/Master-Thesis/Programme/GAMA/feibra_ZEIT_punkte.csv"
	file delivery_time_file <- csv_file("Path/to/your/csv.csv",";",true);
	
	//this file should contain a connected detailed line shapefile for the deliverer to follow on foot.
	//Example path: "D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502601/Straßennetzwerk_Lambert_31287_Gebiet_502601_für den deliverer_von_OSM.shp"
	file street_for_deliverer <- file("Path/to/your/shapefile.shp");
	
	//this file should contain the csv you want to save the results to
	//Example path: "D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502601/Results/Results_Batch_1.csv"
	file my_csv_file <-  csv_file("Path/to/your/csv.csv", ",",float,true);
	
	// setting the geometry of the world to the envelope of the street shapefile
	geometry shape <- envelope(street_shapefile);
	
	// setting the steps to 1s
	float step <- 1 #s;
	
	//creating a graph for the vehicle to follow
 	graph the_graph;
 
 	//creating a graph for the deliverer to follow
 	graph the_walking_graph;	
 	
 	// setting the movement speed of the deliverer outside the car
 	float movespeed <- 4 #km / #h;
 	
 	// setting the default movement speed of the vehicle 
 	float vehiclespeed <- 50 #km / #h;
 	
 	// values for time tracking
 	float total_time <- 0;
 	float total_time_min <-0;
 	float total_time_h <-0; 
 	float time_for_delivery;
 	float time_for_delivery_h;
 	float time_without_delivery_h;
 	float time_for_vehicle_interaction;
 	float time_for_vehicle_interaction_h;
 	float total_time_without_driving_from_or_to_the_area; 
 	float time_to_area<-0.0;
	float time_from_area<-0.0;
 	
 	// setting the time it takes to enter or exit the vehicle and refill the parcels 	
 	float enter_exit_time <- 30#s;
 	
 	// load the time data from the CSV into a matrix
 	matrix time_matrix <- delivery_time_file;
 	
 	// store the length of the matrix
  	int matrix_length <- length(time_matrix);
  	
  	// maximum number of parcels the deliverer can carry
  	int maximum_parcel_carry_capacity <- 2000;
  	
  	// fixed distance at which the deliverer decides to drive instead of walking
  	int distance_until_driving_fixed <-500;
  	
  	// distance at which the deliverer decides to drive instead of walking may change for failsave reasons therefor a fixed and changeable variable exist
 	int distance_until_driving <- distance_until_driving_fixed;
 	
 	// bool to stop each simulation during batch simulation
 	bool fertig <- false;
 	
 	//number that indicates the number of the run in a batch simulation - used to use different adresses as first adress
 	int simulation_counter <-  0 ;
 	
 	//bool for failsave reasons
 	bool endgame <-false;
 	
 	// this map is needed for debugging reasons
	map<point, string> point_to_address_name <- map([]);
	
	//this geometry will be saved as output
	list<geometry>  total_path_traveled;
	
	// this value saves the travelled path 
	float  total_path_traveled_length;
	
	// saves the location of the deliverers previous location for travelled path calculation
	point previous_deliverer_location;

 	//this list is used to recognize and force points of entry to adresses onto the road for the deliverer.
 	list<adress> no_way_found_list <- [];
 	
 	// Define a search radius 
 	float search_radius <- 100.0; 
 	
 	
 	
 	
    init {    	
    	// creating streets from the shapefile with the maximum driving speed as variable maxspeed
    	//Example: [maxspeedT::read('VMAX_CAR_T'),maxspeedB::read('VMAX_CAR_B')]
        create street from: street_shapefile with: [maxspeedT::read('column_of_maximum_speed'),maxspeedB::read('column_of_maximum_speed_other direction')];

        // create walkable streets from the shapefile for the deliverer
        create walkablestreet from: street_for_deliverer;
        
        //creating a graph for the vehicle to follow
        ask street {
        do setup_graph;
    }
    	// creating a weigth map for the street graph using the length of each street
    	map<street,float> weights_map_for_street <- street as_map (each::each.shape.perimeter);
    	
    	//creating the street graph as directed graph with weigths
    	the_graph <- directed (as_edge_graph(street) with_weights weights_map_for_street) ;
    	
    	
        //creating a graph for the deliverer to follow using the length of each street
        map<walkablestreet,float> weights_map_for_walkablestreet <- walkablestreet as_map (each::each.shape.perimeter);       
        the_walking_graph <- as_edge_graph(walkablestreet) with_weights weights_map_for_walkablestreet;
               
        // creating the starting point where the deliverer and vehicle are created
        create startpunkt from: start_point;
        
        //creating buildings --> optional
        create building from: building_shapefile;
        
        // creating adresses from the shapefile with the number of households as variable households
        // Example: with: [households::read('HAUSHALTE')];
        create adress from: adress_shapefile with: [households::read('column_of_households')];
        
        // creating the deliverer with the start at the starting point
        create deliverer number: 1{location <- any (start_point);}
        
        // creating the vehicle with the start at the starting point
        create vehicle number: 1 {location <- any (start_point);}    
        
        //find the closest street to the staring point and position the vehicle on the closest point to the starting point on the street
        ask vehicle {
        	do find_start;        	        	
        }
              
        //creating entry points for the deliverer and parking spots for the vehicle      
        ask adress {
        	if (households = 0){
        		do die;        		
        		}
        		do find_closest_point_on_foot;
        		do find_closest_point_on_driveable_street;        	
 			    point_to_address_name[closest_point_on_driveable_street] <- name;
 			   
        		
        		
        		
 			    
        } 
       // find roads on the walking_graph which contain an entry point and have no path(are not connected) to the deliverer
       ask adress {
       	 path path_to_adress <- path_between(the_walking_graph, any(deliverer), closest_point_on_foot);
			
			if path_to_adress = nil {
			
			no_way_found_list <- no_way_found_list +self;
			ask walkablestreet {
				if distance_to(self.shape,  myself.closest_point_on_foot) < 20 {
					isolate <- true;
					
				}
			}
			
       }
       
       }
       // move the points to a connected segment of the graph
       do adjust_closest_points_to_graph;
       //check if there are still unconnected entry points
               ask adress {
       	 path path_to_adress <- path_between(the_walking_graph, any(deliverer), closest_point_on_foot);
			
			if path_to_adress = nil {
			write self.name + "Error! some Paths cant be found!";
			no_way_found_list <- no_way_found_list +self;
			
       }
       
       }
        write "Number of adresses to be delivered: " + length (adress);
          
        //moves the deliverer to the car
     	ask deliverer {
     		// uncomment the next line if you want to add the time of the deliverer walking to the car to your total time
     		//do gotocar;
     		location <- any (vehicle).location;
     	}    	    	
    }// end of init
    
    // Action to retrieve nearby segments from the graph within a search radius
action get_nearby_segments(point current_point) {
    list<geometry> nearby_segments <- [];
    
    // Loop through all segments in the graph
    ask walkablestreet {
    	if not isolate{
        geometry segment <- self;
        // Check if the distance to the segment is within the search radius
        if (current_point distance_to segment < search_radius) {
            nearby_segments <- nearby_segments + segment;
        }
        }
    }
    
    return nearby_segments;
}


// Action to adjust points to their closest location on the graph
action adjust_closest_points_to_graph {
    // Iterate over all addresses in the no_way_found_list
    loop addr over: no_way_found_list {
    	
        point problematic_point <- addr.closest_point_on_foot;
        float min_distance <- 1e9;
        point closest_point_on_segment <- nil;

        // Retrieve nearby segments for efficiency
        list<geometry> walking_segments <- get_nearby_segments(problematic_point);

        // Iterate over nearby segments
        loop segment over: walking_segments {
            list<point> segment_points <- segment.points;

            // Ensure the segment has at least two points to form a line
            if length(segment_points) > 1 {
                loop i from: 0 to: (length(segment_points) - 2) {
                    point p1 <- segment_points[i];    

                    // Calculate the distance to the projected point
                    float distance <- problematic_point distance_to p1;

                    // Update the minimum distance and closest point if necessary
                    if distance < min_distance {
                        min_distance <- distance;
                        closest_point_on_segment <- p1;
                    }
                }
            }
        }

        // Adjust the closest_point_of_foot to the nearest point of a segment on the graph
        if closest_point_on_segment != nil {
            addr.closest_point_on_foot <- closest_point_on_segment;
        }
    }
}









    //counts the seconds for each tick
    reflex time_passing {
    	total_time <- total_time +1;
    }
    
    //tracks the distance and geometry traveled
    reflex track_deliverer_path {
        // Get the current position of the deliverer
        point current_position <- one_of(deliverer).location;

        // If there's a previous position, create a segment and add it to the total path
        if previous_deliverer_location != nil {
            geometry segment <- line([previous_deliverer_location, current_position]);
            if total_path_traveled = nil {
                total_path_traveled <- segment;
            } else {
                total_path_traveled <- total_path_traveled +segment;
            }
            total_path_traveled_length <- total_path_traveled_length + segment.perimeter;
        }

        // Update the previous position
        previous_deliverer_location <- current_position;
    }
    
    //condition to end the simulation after all adresses are deleted 
    reflex end_simulation when: length(adress) = 0 { 
    	
    	ask deliverer {
    		if on_wayout = false and incar =false{
    			
    		do return_to_vehicle; 
    		 		
    		}
    		
    	if incar {
    	ask vehicle {
    		//drives the vehicle to the starting position
    		do return_to_base;
    		if location = closest_starting_point_on_street.location {
			do die;
		}
		
		}
		
		}
    	}
    	if length(vehicle) = 0 {
    		//stops the simulation after the vehicle was killed and saves the data
    		total_time_min <- total_time / 60;
    		total_time_h <- total_time / 60 / 60;
    		time_for_delivery_h <- time_for_delivery / 60 / 60;
    		time_without_delivery_h <- total_time_h -time_for_delivery_h;
    		time_for_vehicle_interaction_h <- time_for_vehicle_interaction  / 60 / 60;
    		total_time_without_driving_from_or_to_the_area <- (total_time - time_to_area) - time_from_area;
    		write "Time in seconds: " + total_time with_precision 2; 
    		write "Time in minutes: "+total_time_min with_precision 2; 
    		write "Time in hours: "+total_time_h with_precision 2;
    		write "Time in seconds from time table: "+time_for_delivery with_precision 2;
    		write "Time in hours from time table: "+time_for_delivery_h with_precision 2;
    		write "Time in hours without time table: "+time_without_delivery_h with_precision 2;
    		write "Time for the way to the area in seconds: " + time_to_area;
    		write "Time for the way out of the area in seconds: " + time_from_area;
    		write "Time for vehicle interaction: "+time_for_vehicle_interaction_h with_precision 2;
    		write "Time inside the area: " +total_time_without_driving_from_or_to_the_area with_precision 2;
    		save [simulation_counter,total_path_traveled_length, distance_until_driving_fixed,time_for_vehicle_interaction_h, time_for_delivery_h,total_time, time_to_area, time_from_area,total_time_without_driving_from_or_to_the_area] type: csv to: my_csv_file rewrite: false;  
    		
    		// save the path as a shapefile at the end of the simulation    
    		//Example path: "D:/Master/Master-Thesis/Programme/GAMA/Daten für Gebiet 502601/Path.shp"
            if (total_path_traveled != nil) {           	
            save total_path_traveled to: "Path/to/your/shapefile.shp" type: "shp";
            write "Path saved as shapefile: shapefile.shp";
            write "Total path length: " + total_path_traveled_length + " meters.";
        	}    	   	
    		fertig <- true;
    	}
		}
    
}

species building {
	//species buildings serves only for visual reasons
	rgb color <- #gray;
    aspect base {
        draw shape color: color;
    }
}


species walkablestreet {
	//species walkablestreet are all streets that can be accessed by foot
	 rgb color <- #purple ;
	 point first_vertex <- nil;
     point last_vertex <- nil;
     graph street_graph;
     bool isolate <-false;
	 
     aspect base {
   		draw shape color: color ;
   	}
   	
}



species street {
   //species street contains only the streets a car can drive on
   float maxspeedT;
   float maxspeedB;
   bool in_endgame <- false;
   rgb color <- #blue ;
   float street_width <- 2;
   bool arrow_green <- false;
   bool arrow_red <- false;
 
  
   
   aspect base {   	
   	draw shape color: color width: street_width;   	
    if (maxspeedT = -1) {
        do draw_arrow_along(shape.points, #red);
    } else if (maxspeedB = -1) {
        do draw_arrow_along(shape.points, #green);
    } 
   }

   //this action creates streets, so the species can be used as directed graph and one-way streets are considered   
action setup_graph {
    // Check if the road is a one-way street
    if (maxspeedT = -1) {
        // If only backward direction is allowed (one-way street)
        color <- #red;
        
        // Reverse the geometry of the road to match the backward direction
        shape <- polyline(reverse(shape.points));
        
    } else if (maxspeedB = -1) {
        // If only forward direction is allowed (one-way street)
        color <- #green;
        
        // No need to reverse geometry, as it aligns with the forward direction
    } else if (maxspeedT > -1 and maxspeedB > -1) {
        // If the street is bidirectional, create a duplicate street in reverse direction
        create street {
            shape <- polyline(reverse(myself.shape.points));
            color <- myself.color;        
          	}
          	
          	}
            
        }
        
        
  //this action visualises the direction of one way streets          
 action draw_arrow_along(list<point> points, rgb color) {
    float arrow_spacing <- 5.0;  
    float arrow_size <- 3;     
    float wing_offset <- 3;    
    
    // Traverse the points to place arrows at intervals
    float accumulated_distance <- 0.0;
    loop i from: 0 to: (length(points) - 2) {
        point start <- points[i];
        point end <- points[i + 1];
        
        // Calculate the distance between points
        float segment_length <- distance_to(start, end);
        accumulated_distance <- accumulated_distance + segment_length;
        
        // Draw an arrow if we've reached the arrow spacing threshold
        if accumulated_distance >= arrow_spacing {
            accumulated_distance <- 0.0;  // Reset distance accumulator
            
            // Draw main line segment of the arrow
            //draw line([start, end]) color: color;
            
            // Calculate the direction vector
            point direction <- normalize_vector(start, end);
            
            // Calculate the arrowhead points
            point arrow_tip <- [end.x - 0.5 * direction.x, end.y - 0.5 * direction.y];
           point left_wing <- [arrow_tip.x - (wing_offset + arrow_size) * direction.x + arrow_size * (-direction.y), 
                                arrow_tip.y - (wing_offset + arrow_size) * direction.y + arrow_size * direction.x];
            point right_wing <- [arrow_tip.x - (wing_offset + arrow_size) * direction.x - arrow_size * (-direction.y), 
                                 arrow_tip.y - (wing_offset + arrow_size) * direction.y - arrow_size * direction.x];
            // Draw the arrowhead wings
            draw line([arrow_tip, left_wing]) color: color;
            draw line([arrow_tip, right_wing]) color: color;
        }
    }
}

// Helper function to normalize a vector between two points
action normalize_vector(point p1, point p2) {
    float dx <- p2.x - p1.x;
    float dy <- p2.y - p1.y;
    float length <- sqrt(dx^2 + dy^2);
    return [dx / length, dy / length];
}   
   
}

species adress {
   // this species contains the adresses
   int households;
   int recparcels;
   rgb color <- #red;
   float min_distance_to_adress <- 1e9;
   point closest_street_location;
   point first_vertex;
   point last_vertex;
   graph street_graph;
   point temp_closest_point_on_foot;
   point temp_closest_point_on_driveable_street;
   point closest_point_on_driveable_street;
   point closest_point_on_foot;
   float min_point_distance <- 1e9;
   float min_driveable_point_distance <-1e9;
   bool considered_on_way_to_vehicle <- false; 
   list<street> streets_for_endgame;
   point closest_point_on_walkable_path;
   point closest_point_on_foot_to_point_on_driveable_street;

   aspect base {
        draw circle(2) color: #red;
        draw circle(2) at: closest_point_on_foot color: #green; 
        draw circle(2) at: temp_closest_point_on_driveable_street color: #purple; 
        
        
   } 

//action to find the point on the closest road(for the deliverer) which ist closest to the adress - resembles a garden gate
   action find_closest_point_on_foot {
   	    walkablestreet nearest_street <- nil;
        float min_street_distance <- 1e9;
        list<point> points_to_choose <-nil; 
        list<walkablestreet> excluded_streets <-nil;       
        //repeated 10 times - the 10 closest streets are considered              
        loop times: 10{              
       ask walkablestreet {
       	   if (not (self in excluded_streets)){
       		
            	float distance_to_adress <- self.location distance_to myself.location;  
            	if (distance_to_adress < min_street_distance) {
                	min_street_distance <- distance_to_adress;
                	nearest_street <- self;
                geometry street_geometry <- self.shape;
                list<point> street_geometry_points_list <- street_geometry.points;
                myself.first_vertex <- street_geometry.points[0].location;
                myself.last_vertex <- street_geometry.points[length(street_geometry_points_list) - 1].location;
                myself.street_graph <- as_edge_graph(nearest_street);
                }
            
            }
        }
        
        if (nearest_street != nil) {
     	  excluded_streets <- excluded_streets + nearest_street;
     	    } 	

     	//create a not visible walker which creates points along a road (1m steps) - the closest point to the adress will be saved 
     	create walker number: 1 {
     		location <- myself.first_vertex.location;
     		}
     		
     	ask walker {
     		bool stillgoing <- true;
     		loop while: stillgoing {
     			do goto target:myself.last_vertex on: myself.street_graph speed: 1.0;
     			float distance_adress <- location distance_to myself.location;
     			
     			if (distance_adress < myself.min_point_distance) {
            	    myself.min_point_distance <- distance_adress;
            	    myself.temp_closest_point_on_foot <- self.location;
           	     }
     			if location = myself.last_vertex.location {
     				stillgoing <- false;
     				points_to_choose <- points_to_choose + myself.temp_closest_point_on_foot;
     				}
     			}
     		do die;
     		}
     	      	       	if self.name = "adress283"{

     	}
     	
     	   min_point_distance <- 1e9;	
     	   min_street_distance <- 1e9;
     	   
 
     	}
     	
     	
     	
	//choosing the closest of all 10 closest points
	closest_point_on_foot <- with_min_of(points_to_choose, each distance_to location);

	}
	
//action to find the point on the closest street (for the vehicle) which ist closest to the adress - resembles a parking spot
   action find_closest_point_on_driveable_street {
   	    street nearest_street <- nil;
        float min_street_distance <- 1e9;
        list<point> points_to_choose <-nil; 
        list<street> excluded_streets <-nil;       
        //repeated 10 times - the 10 closest streets are considered              
        loop times: 10{              
       ask street {
       	   if (not (self in excluded_streets)){
       		
            	float distance_to_adress <- self.location distance_to myself.location;  
            	if (distance_to_adress < min_street_distance) {
                	min_street_distance <- distance_to_adress;
                	nearest_street <- self;
                geometry street_geometry <- self.shape;
                list<point> street_geometry_points_list <- street_geometry.points;
                myself.first_vertex <- street_geometry.points[0].location;
                myself.last_vertex <- street_geometry.points[length(street_geometry_points_list) - 1].location;
                myself.street_graph <- as_edge_graph(nearest_street);
                }
            
            }
        }
        
        if (nearest_street != nil) {
     	  excluded_streets <- excluded_streets + nearest_street;
     	    } 	
     	//create a not visible walker which creates points along a street (1m steps) - the closest point to the adress will be saved 
     	create walker number: 1 {
     		location <- myself.first_vertex.location;
     		}
     		
     	ask walker {
     		bool stillgoing <- true;
     		loop while: stillgoing {
     			do goto target:myself.last_vertex on: myself.street_graph speed: 1.0;
     			float distance_adress <- location distance_to myself.location;
     			
     			if (distance_adress < myself.min_driveable_point_distance) {
            	    myself.min_driveable_point_distance <- distance_adress;
            	    myself.temp_closest_point_on_driveable_street <- self.location;
           	     }
     			if location = myself.last_vertex.location {
     				stillgoing <- false;
     				points_to_choose <- points_to_choose + myself.temp_closest_point_on_driveable_street;
     				}
     			}
     		do die;
     		}
     	 
     	
     	   min_point_distance <- 1e9;	
     	   min_street_distance <- 1e9;
     	   
     	   
     	}
     	
     	
     	
	//choosing the closest of all 10 closest points
	closest_point_on_driveable_street <- with_min_of(points_to_choose, each distance_to location);

	}	
	
}

species walker skills:[moving] {
            //this species serves only the purpose of finding closest points on lines
    }


species deliverer skills:[moving] {
	//this species is only one agent and contains the deliverer
	float dspeed <- movespeed;
	int maxparcels <- maximum_parcel_carry_capacity;
	int parcels <- maxparcels;
	bool incar <- false;
	rgb color <- #black;
	point the_target <- nil ;
	point the_target_entrance <- nil;
	string the_target_name <- nil;
	float distance_to_next_adress <- nil;
	int needed_parcels;
	bool in_entrance <- false;
	point exitpoint;
	bool on_wayout <- false;
	bool loopstopper <- false;
	path show_me_the_path;
	list<string> already_considered <- [];
	bool after_delivery <- false;
	
	
	aspect base{
				// loads a svg file to visualise the deliverer 
		// Example path:"D:/Master/Master-Thesis/Programme/GAMA/deliverer.svg"
		draw file("D:/Master/Master-Thesis/Programme/GAMA/deliverer.svg") at: {location.x, location.y, 10} size: 10 color: color rotate:180;
			//use the following line without "//" to show the next targets entry point
		//draw circle(100) at:the_target_entrance color: #orange; 
 
	} 
	//when in a vehicle the deliverer drives 
	reflex drive_on when: incar {
		ask vehicle {
			if length(adress) !=0{
			do drive;			
			}
		}
	}
	
	 // initial reflex to find new target
    reflex initial_target when: the_target = nil and not incar and not on_wayout{
        do find_nearest_adress;
    }
	// Method to set the nearest address(actually the nearest entry point) as destination
	action find_nearest_adress {
		
		if length(adress) != 0 and the_target = nil{
	
        adress nearest_adress <- nil;        
        // setting a very big number as min_distance
        float min_distance <- 1e9;  
        
        ask adress {
           path path_to_adress <- path_between(the_walking_graph, myself.location, closest_point_on_foot.location);
			
			if path_to_adress != nil {
				
			    float total_length <- 0.0;
 			   list<geometry> segments <- path_to_adress.segments;

 			   loop segment over: segments {
 			       float segment_length <- segment.perimeter;
 			       total_length <- total_length + segment_length;
  			  }

  			  float distance_to_adress <- total_length;

 			   if (distance_to_adress <= min_distance) {
 			       min_distance <- distance_to_adress;
 			       nearest_adress <- self;
  			      myself.needed_parcels <- self.households;
 			   }
			}
			
			
			}
			
		
        if (nearest_adress != nil) {
            the_target <- nearest_adress.location; 
            the_target_name <- nearest_adress.name;
            the_target_entrance <- nearest_adress.closest_point_on_foot.location; 
			
        
        
        
        
        show_me_the_path <- path_between(the_walking_graph, location, nearest_adress.closest_point_on_foot.location);
        float distance_to_next_adress_on_path <- 0;
		list<geometry> segments <- show_me_the_path.segments;
		loop segment over: segments {
            float segment_length <- segment.perimeter;  
            distance_to_next_adress_on_path <- distance_to_next_adress_on_path + segment_length; 
        }
        distance_to_next_adress <- distance_to_next_adress_on_path;
        
        if endgame {       	
        	distance_until_driving <-distance_to_next_adress +1;
       
        }
}else{
	write "Error: no nearest adress!";
	do return_to_vehicle;
}
        
        
        }
	}
	//this action finds the nearest adress an the way to the car and prevents a loop triggered by changing targets -->loopstopper
	action find_nearest_adress_during_return {
		if length(adress) != 0 and the_target = nil{
        adress nearest_adress <- nil;
        float min_distance <- 1e9;  
        
        ask adress {
           path path_to_adress <- path_between(the_walking_graph, myself.location, closest_point_on_foot);
			
			if path_to_adress != nil {
			    float total_length <- 0.0;
 			   list<geometry> segments <- path_to_adress.segments;

 			   loop segment over: segments {
 			       float segment_length <- segment.perimeter;
 			       total_length <- total_length + segment_length;
  			  }

  			  float distance_to_adress <- total_length;

 			   if (distance_to_adress <= min_distance) {
 			       min_distance <- distance_to_adress;
 			       nearest_adress <- self;
  			      myself.needed_parcels <- self.households;
 			   }
			}
			
			
			}
        if (nearest_adress != nil) {
            the_target <- nearest_adress.location; 
            the_target_name <- nearest_adress.name;
            the_target_entrance <- nearest_adress.closest_point_on_foot.location; 
        
        
 		path path_to_adress <- path_between(the_walking_graph, location, nearest_adress.closest_point_on_foot.location);
           float distance_to_next_adress_on_path <- 0;
		list<geometry> segments <- path_to_adress.segments;
		loop segment over: segments {
            float segment_length <- segment.perimeter;  // Berechne die Länge des Segments
            distance_to_next_adress_on_path <- distance_to_next_adress_on_path + segment_length;  // Füge die Segmentlänge zur Gesamtlänge hinzu
        }
        distance_to_next_adress <- distance_to_next_adress_on_path;
        if endgame {
        	if not after_delivery{
        	distance_until_driving <-distance_to_next_adress +1;
        	
        	}
        }
        if needed_parcels <= parcels {
        	loopstopper <-true;
        }
        }else{
        	write "Error: no nearest adress during return!";
        	}
	}
	}
	
	
	//considering return to vehicle
	//return to vehicle if the deliverer has not enough parcels or the next adress is too far away - the deliverer can not return to the vehicle until he is back at the enter point 
	reflex considering when: incar = false {
		        	
		if ((needed_parcels > parcels or distance_to_next_adress > distance_until_driving) and on_wayout=false){
			the_target <- nil; 
            the_target_name <- nil;
            the_target_entrance <- nil; 
			do return_to_vehicle;
		}
		else {
			do walk;
		}
		
	}
	
	//action to return the deliverer to the vehicle	
	action return_to_vehicle {
        float distance_to_vehicle <- location distance_to any (vehicle).vehicle_substitute.location;
        point vehicle_location <- location of one_of (vehicle).vehicle_substitute;
        point true_vehicle_location <- location of one_of (vehicle);
       
 if (distance_to_vehicle > 0) {
 	
        do goto target: any(vehicle).vehicle_substitute.location on: the_walking_graph speed: dspeed;
        
    //there is a bug where the deliverer does not move after entering the walking_graph after the last adress had been killed - this if forces a second move and this fixes the problem 
    if length(adress) = 0{    
    do goto target: any(vehicle).vehicle_substitute.location on: the_walking_graph speed: dspeed; 
    total_time <- total_time +1;
    }   
    //}

    
}

        	//find a new target as there might be an opportunity on the way to the vehicle
        	if length (adress) !=0 {
        	do find_nearest_adress_during_return;
        	
        	}

        	// checking if deliverer is at vehicle and refill the parcels
        	if (location = vehicle_location) {
            	parcels <- maxparcels;
            	
            	do find_nearest_adress;
            	if length(adress) = 0{
            		location <- true_vehicle_location;
            		incar<-true;
            	}
        	}
        	//if the deliverer enters the car time will pass, deliverer enters if next adress is too far away to walk
        	if (distance_to_next_adress >= distance_until_driving and location = vehicle_location){
        		incar <-true;
        		location <- true_vehicle_location;
        		loopstopper <- false;  
        		total_time <- total_time + enter_exit_time;
        		time_for_vehicle_interaction <- time_for_vehicle_interaction+ enter_exit_time;       		 						
        	}
        	
	}
	 // movement to the target on the walkable streets
    action walk {
    	if (distance_to_next_adress < distance_until_driving and parcels >= needed_parcels and  in_entrance = false) {
    		
    		do goto target:the_target_entrance on: the_walking_graph speed:dspeed;
    		
        
        // checking if deliverer is at target location and entering the entrance
        if location = the_target_entrance {
        	        	in_entrance <- true;
        	        	
        }
        }
        //moving without a graph in a straigth line
        if in_entrance = true and on_wayout = false{
        	do goto target: the_target speed: dspeed;
        	
        }
        //returning to the entrance after delivery
        if in_entrance = true and on_wayout = true {      	
        	do goto target: exitpoint speed: dspeed;
        	//resetting 
        	if location = exitpoint {
        		exitpoint <- nil;
        		on_wayout <- false;
        		in_entrance <- false;
        		after_delivery <- true;
 		       	after_delivery <- false;
 		       	the_target <- nil;
   		     	the_target_name <- nil;
        		do find_nearest_adress;
        	}
        }
        if (location = the_target) {
        	// deliver parcels if at adress
        	do deliver_parcel;  
        }
    }
    
    action walk_without_checking {
    	//should be used when going to a adress without searching for better (closer adresses) to prevent endles loops
    	if (distance_to_next_adress < distance_until_driving and parcels >= needed_parcels and  in_entrance = false) {
    		
    		do goto target:the_target_entrance on: the_walking_graph speed:dspeed;
            if location = the_target_entrance {
        	        	in_entrance <- true;
        	        	
        }
        }
        if in_entrance = true and on_wayout = false{
        	do goto target: the_target speed: dspeed;
        	
        }
        if in_entrance = true and on_wayout = true {
        	do goto target: exitpoint speed: dspeed;
        	if location = exitpoint {
        		exitpoint <- nil;
        		on_wayout <- false;
        		in_entrance <- false;
        		after_delivery <- true;
 		       	after_delivery <- false;
 		       	the_target <- nil;
   		     	the_target_name <- nil;
    	        if endgame {        	        	
	        	endgame <-false;
       			}
        		}
        		
        }
        if (location = the_target) {
            do deliver_parcel;  
        }
    }
    
    // deliver parcels and delete adress 
    action deliver_parcel {
    	// set exitpoint so the deliverer can return to the street
    	exitpoint <- the_target_entrance;
    	on_wayout <-true;
        ask adress {
        	if (self.location = myself.location){
            myself.parcels <- myself.parcels - households;
            recparcels <- recparcels + households;
            //looking inside the time matrix on position of the number of households (-1 because 1 household ist at line 0 )
            string stringfrommatrix <- time_matrix[1,households-1];
            float floatfrommatrix <- float(stringfrommatrix);
            //time is in min --> conversion to s
            float time_per_address <- floatfrommatrix * 60;          
            time_for_delivery <- time_for_delivery + time_per_address;
            //adding the time
            total_time <- total_time + time_per_address;
			// deleting of adress after delivery
            if (recparcels >= households) {
                do die;
            }
            
            }
        }
        
        // find next adress
        the_target <- nil;
        the_target_name <- nil;
        if endgame {        	        	
        	endgame <-false;
        }
       }
    //gotocar serves the purpose of positioning the deliverer at the car on the street. the time walking there will be added.
   action gotocar {
   	bool ongoing <- true;
   	loop while: ongoing {
   		float move_distance <- min(location distance_to any(vehicle), dspeed * step);
   	    point carlocation <- any (vehicle);
     point move_step <- (location + ((carlocation - location)*(move_distance / (location distance_to carlocation))));
       	        
        //moving the agent to the next location
        location <- move_step;
        total_time <- total_time +1;
        if location = carlocation {
        	ongoing <- false;
        	incar <-true;
        	total_time <- total_time + enter_exit_time;
        	time_for_vehicle_interaction <- time_for_vehicle_interaction+ enter_exit_time;
        }
         
   }
   
   }
   
  
   
}






species vehicle skills: [moving] {
	//this species is only one agent and contains the vehicle 
	float vspeed <- vehiclespeed;
	rgb color <- #orange;
	point next_parking_spot <- nil;
	float vehicle_Start_min_dist <- 1e9;
	float vehicle_parking_spot_min_dist <- 1e9;
	point closest_street_location;
    point first_vertex_vehicle;
    point last_vertex_vehicle;
    graph street_graph_vehicle; 
    float min_start_point_distance <- 1e9;
	point temp_closest_start_point_on_street;
	point closest_starting_point_on_street;
	list<point> possible_parking_spots;
	list<geometry> exact_possible_parking_spots;
	list<point> choosen_fallback_parking_spots;
	bool is_forward <- true;
	point previous_location <- nil;
	bool first_parking_spot_reached <- false;
	list<street>streets_for_endgame;
	path path1;
	path path2;
	string next_parking_spot_name;
    point vehicle_substitute <- nil;  
    
	
	
	aspect base{
		// loads a svg file to visualise the vehicle 
		// Example path:"D:/Master/Master-Thesis/Programme/GAMA/car.svg"
		draw file("path/to/your/svg.svg") size: 10 color: color rotate:180; 
		if (next_parking_spot != nil) {
            draw circle(10) at: {next_parking_spot.x,next_parking_spot.y, 9} color: #blue border: #black;  // Darstellung des nächsten Parkplatzes
            loop i over:choosen_fallback_parking_spots{
            draw circle (5) at: i color: #green;
            
            }
        draw path1.shape color: #green width: 5;
        draw path2.shape color: #purple width: 5;
        }
	} 
	
	
	//finding the nearest road to the starting point and positioning the vehicle there
	action find_start {
		street nearest_starting_street <- nil;
        float min_starting_street_distance <- 1e9;
        list<point> points_to_choose_start; 
        list<street> excluded_starting_streets;       
        //looping 4 times to consider the 4 nearest roads               
        loop times: 4{              
       ask street {
       	   if (not (self in excluded_starting_streets)){
       		
            	float distance_to_start <- self.location distance_to any(start_point);  
            	if (distance_to_start < min_starting_street_distance) {
                	min_starting_street_distance <- distance_to_start;
                	nearest_starting_street <- self;
                geometry street_geometry_vehicle <- self.shape;
                list<point> street_geometry_vehicle_points_list <- street_geometry_vehicle.points;
                myself.first_vertex_vehicle <- street_geometry_vehicle.points[0].location;
                myself.last_vertex_vehicle <- street_geometry_vehicle.points[length(street_geometry_vehicle_points_list) - 1].location;
                myself.street_graph_vehicle <- as_edge_graph(nearest_starting_street);
                }
            
            }
        }
        
        if (nearest_starting_street != nil) {
     	  excluded_starting_streets <- excluded_starting_streets + nearest_starting_street;
     	    } 	
     	 
     	create walker number: 1 {
     		location <- myself.first_vertex_vehicle.location;
     		}
     		
     	ask walker {
     		bool stillgoing <- true;
     		loop while: stillgoing {
     			do goto target:myself.last_vertex_vehicle on: myself.street_graph_vehicle speed: 1.0;
     			float distance_start <- location distance_to any(start_point);
     			
     			if (distance_start < myself.min_start_point_distance) {
            	    myself.min_start_point_distance <- distance_start;
            	    myself.temp_closest_start_point_on_street <- self.location;
           	     }
     			if location = myself.last_vertex_vehicle.location {
     				stillgoing <- false;
     				points_to_choose_start <- points_to_choose_start + myself.temp_closest_start_point_on_street;
     				}
     			}
     		do die;
     		}
     	 
     	
     	   min_start_point_distance <- 1e9;	
     	   min_starting_street_distance <- 1e9;
     	        	   
     	}
     	    	
     	
closest_starting_point_on_street <- with_min_of(points_to_choose_start, each distance_to any(start_point));
location <- closest_starting_point_on_street.location;
}



// Action to find parking spots
action find_parking_spots {
    possible_parking_spots <- nil;
    exact_possible_parking_spots <- nil;
    vehicle_parking_spot_min_dist <- 1e9;

//     Start at a random parking spot if first parking spot has not been reached
    if (first_parking_spot_reached = false) {    
        next_parking_spot <- adress[simulation_counter].closest_point_on_driveable_street;   
        
        path path_to_next_spot  <- path_between(the_graph, location, next_parking_spot);
     	
     	if (path_to_next_spot = nil) {
            write "No valid path to the first parking spot. Nearest pointwill be chosen.";
            // Find the nearest connected non-oneway street
            street nearest_connected_street <- nil;
            float min_distance <- 1e9;
			list<street> nearest_streets_ini <- [];
			
			ask adress[simulation_counter] {               
                float min_distance;
                loop times: 4 {
                    float nearest_distance <- 1e9;
                    street nearest_street <- nil;
                    ask street {
                        float distance_to_adress <- location distance_to myself.closest_point_on_driveable_street;
                        if (distance_to_adress < nearest_distance and not(self in nearest_streets_ini)) {
                            nearest_distance <- distance_to_adress;
                            nearest_street <- self;
                        }
                    }
                    if nearest_street != nil {
                        nearest_streets_ini <- nearest_streets_ini + nearest_street;
                    }
                
               }
                           
                }   
			
			float end_distance <- 1e9; // Initialize with a large value
            // Iterate over all streets to find a suitable connected endpoint
            loop i over: nearest_streets_ini {
                if (i.maxspeedT > -1 and i.maxspeedB > -1) { // Exclude one-way streets in the wrong direction
                    list<point> street_points <- i.shape.points;

                    // Check the start and end points of the street
                    point start_point_normal_street <- street_points[0];
                    point end_point_normal_street <- street_points[length(street_points) - 1];

                    // Determine the connection distance
 					// Calculate the path distance to the start point
					path path_to_start <- path_between(the_walking_graph, next_parking_spot, start_point_normal_street);
                    float start_distance <- 1e9; // Initialize with a large value
                    if (path_to_start != nil) {
                        start_distance <- 0.0;
                        list<geometry> segments <- path_to_start.segments;
                        loop segment over: segments {
                            start_distance <- start_distance + segment.perimeter;
                        }
                    }

                    // Calculate the path distance to the end point
                    path path_to_end <- path_between(the_walking_graph, next_parking_spot, end_point_normal_street);
                    
                    if (path_to_end != nil) {
                        end_distance <- 0.0;
                        list<geometry> segments <- path_to_end.segments;
                        loop segment over: segments {
                            end_distance <- end_distance + segment.perimeter;
                        }
                    }
                    
                    if (start_distance < min_distance) {
                        min_distance <- start_distance;
                        nearest_connected_street <- self;
                        next_parking_spot <- start_point_normal_street;
                    }

                    if (end_distance < min_distance) {
                        min_distance <- end_distance;
                        nearest_connected_street <- self;
                        next_parking_spot <- end_point_normal_street;
                    }
                }
            }
        }
        
     	         
		if (next_parking_spot in point_to_address_name) {
        next_parking_spot_name <- point_to_address_name[next_parking_spot];
  		  } else {
        next_parking_spot_name <- "Unknown"; 
        
    }   
    
      
    } else {
        // Find possible parking spots near addresses within walking distance
       ask adress {
            path distance_to_enter <- path_between(the_walking_graph, self.closest_point_on_driveable_street, self.closest_point_on_foot);
            float total_length <- 0.0;
            if distance_to_enter != nil{
            list<geometry> segments <- distance_to_enter.segments;
            loop segment over: segments {
                total_length <- total_length + segment.perimeter;
            }

            if total_length < distance_until_driving {
                myself.possible_parking_spots <- myself.possible_parking_spots + closest_point_on_driveable_street;                        
            }
            }
         
        }

        // Calculate the closest parking spot based on path length
        if possible_parking_spots != [] {
            loop i over: possible_parking_spots {
                path path_to_next_spot <- path_between(the_graph, location, i);
if (path_to_next_spot = nil) {
	possible_parking_spots <- possible_parking_spots - i;
}
                if (path_to_next_spot != nil) {
                    float total_length <- 0.0;
                    list<geometry> segments <- path_to_next_spot.segments;

                    // Add the length of all segments except the last one
                    loop segment over: segments {
                        if segment != segments[length(segments) - 1] {
                            total_length <- total_length + segment.perimeter;
                        }
                    }

                    // Handle the last segment only up to the target point
			        if length(segments) > 0 {
  			          geometry last_segment <- segments[length(segments) - 1];
  			          list<point> last_segment_points <- last_segment.points;
   			         point start_node_of_last_segment <- last_segment_points[0];
  			          float segment_distance_to_target <- start_node_of_last_segment distance_to i;

    			        // Add the distance from the start node of the last segment to the target point
    			        total_length <- total_length + segment_distance_to_target;
    			    }

                    // Update the next parking spot if this one is closer
                    if (total_length < vehicle_parking_spot_min_dist) {
                        vehicle_parking_spot_min_dist <- total_length;
                        next_parking_spot <- i;
                        next_parking_spot_name <- point_to_address_name[i];
                    }
                }
            }
           
        }
 
		if possible_parking_spots = [] {
		
			endgame <-true;
            // If no possible parking spots found, look for fallback spots
            list<point> fallback_parking_spots <- nil;
            list<point> fallback_parking_spots_before <- nil;
            // Find the 4 nearest streets to each closest_point_on_foot
            ask adress {
                list<street> nearest_streets <- [];
                float min_distance;
                loop times: 4 {
                    float nearest_distance <- 1e9;
                    street nearest_street <- nil;
                    ask street {
                        float distance_to_adress <- location distance_to myself.closest_point_on_foot;
                        if (distance_to_adress < nearest_distance and not(self in nearest_streets)) {
                            nearest_distance <- distance_to_adress;
                            nearest_street <- self;
                        }
                    }
                    if nearest_street != nil {
                        nearest_streets <- nearest_streets + nearest_street;
                    }
                
                           
                }     
                // Find the closest points on these streets to the closest_point_on_foot
                loop each_street over: nearest_streets {
                    float min_fallback_distance <- 1e9;
                    point best_fallback_spot <- nil;
                    list<point> street_points <- each_street.shape.points;
                    
                   loop street_point over: street_points {
                        path path_to_foot <- path_between(the_walking_graph, street_point, closest_point_on_foot);
                        if path_to_foot != nil {
                            float path_length <- 0.0;
                            list<geometry> segments <- path_to_foot.segments;
                            loop segment over: segments {
                                path_length <- path_length + segment.perimeter;
                            }
                             if path_length < min_fallback_distance {
                                min_fallback_distance <- path_length;
                                best_fallback_spot <- street_point;
                            }                                                       
                            }
                            }
                            if best_fallback_spot != nil {
                        fallback_parking_spots_before <- fallback_parking_spots_before + best_fallback_spot;
                    }
                            }
                            
                   loop each_pos_point over: fallback_parking_spots_before {
                   	point best_fallback_spot <- nil;
                    float min_fallback_distance <- 1e9;
                     path path_to_foot <- path_between(the_walking_graph, each_pos_point, closest_point_on_foot);
                        if path_to_foot != nil {
                            float path_length <- 0.0;
                            list<geometry> segments <- path_to_foot.segments;
                            loop segment over: segments {
                                path_length <- path_length + segment.perimeter;
                            }
                             if path_length < min_fallback_distance {
                                min_fallback_distance <- path_length;
                                best_fallback_spot <- each_pos_point;
                            }                                                       
                            }                         
                            if best_fallback_spot != nil {
                        fallback_parking_spots <- fallback_parking_spots + best_fallback_spot;
                    }
                    }
                    }
                            
                    // Select the fallback parking spot that is closest to the vehicle's current location
            float min_distance_to_vehicle <- 1e9;
            point chosen_fallback_spot <- nil;
                    
                    loop fallback_spot over: fallback_parking_spots {
                    path path_to_vehicle <- path_between(the_graph, location, fallback_spot);
                    
                     if path_to_vehicle != nil {
                    float total_path_length <- 0.0;
                    list<geometry> segments <- path_to_vehicle.segments;
					
					if segments != []{	
                    // Calculate the length of the path
                    loop segment over: segments {
                    	if segment != segments[length(segments) - 1] {
                        total_path_length <- total_path_length + segment.perimeter;
                        }
                    }
                    
                     // Now handle the last segment only up to the target point
                    geometry last_segment <- segments[length(segments) - 1];
                    list<point> last_segment_points <- last_segment.points;
                    point start_node_of_last_segment <- last_segment_points[0];
                    float segment_distance_to_target <- start_node_of_last_segment distance_to fallback_spot;
                    
                     // Add the distance from the start node of the last segment to the target point
                    total_path_length <- total_path_length + segment_distance_to_target;

                    // Check if this path is the shortest path to a fallback spot
                    if total_path_length < min_distance_to_vehicle {
                        min_distance_to_vehicle <- total_path_length;
                        chosen_fallback_spot <- fallback_spot;
                    }
                    }
                }
                    }
                    
              if chosen_fallback_spot != nil {
              	
                next_parking_spot <- chosen_fallback_spot;
                if (next_parking_spot in point_to_address_name) {
        			next_parking_spot_name <- point_to_address_name[next_parking_spot];
  		  			} else {
       				next_parking_spot_name <- "Unknown3";       				
       				}
            } else {
                // If no fallback spot found, choose a random address's closest point as a last resort
                
                    next_parking_spot <- any(adress).closest_point_on_driveable_street;
                    if (next_parking_spot in point_to_address_name) {
        			next_parking_spot_name <- point_to_address_name[next_parking_spot];
  		  			} else {
       				next_parking_spot_name <- "Unknown4";       				
       				}
                
            }
                    }
                                  
                 }  
                 }                                         	
     
  //this action makes the vehicle move, teleports the deliverer to the vehicle and if reached its target creates a substitute point on the nearest road for the deliverer as "car position"
	action drive  {
		float nearest_maxspeed <- 0 #km / #h;
		float min_distance <- 1e9;
        street nearest_street <- nil;
        if next_parking_spot = nil{
        	do find_parking_spots;
          	}
        	
        // find nearest street to set the speed of the vehicle
        ask street {
            float distance_to_vehicle <- self.location distance_to (location);
            if (distance_to_vehicle < min_distance) {
                min_distance <- distance_to_vehicle;
                nearest_street <- self;
            }
        }

        if (nearest_street != nil) {
            nearest_maxspeed <- nearest_street.maxspeedT;
        }
		 vspeed <- nearest_maxspeed #km / #h;
		 
		do goto(target: next_parking_spot, on:the_graph, speed: vspeed);
	ask deliverer {
		location <- myself.location;
	}
		if (next_parking_spot = self.location ) {
			if first_parking_spot_reached =false{
				time_to_area <- total_time;
				
			}
			first_parking_spot_reached <- true;
			next_parking_spot <- nil ;
			// Create a substitute point on the walking graph near the vehicle's location
        	list<point> Possible_points;
        	list<walkablestreet> excluded_walkable_streets ;
        	float min_walkable_street_distance <- 1e9;
        	walkablestreet nearest_walkable_street;
       loop times:4{ 	
        ask walkablestreet {
       	   if (not (self in excluded_walkable_streets)){
       		
            	float distance_to_vehicle <- self.location distance_to myself.location;  
            	if (distance_to_vehicle < min_walkable_street_distance) {
                	min_walkable_street_distance <- distance_to_vehicle;
                	nearest_walkable_street <- self;
                }
            
            }              
        }
        
        if (nearest_walkable_street != nil) {
     	  excluded_walkable_streets <- excluded_walkable_streets + nearest_walkable_street;
     	    }
     	    
     	    }
     	loop i over: excluded_walkable_streets {
     		 	
     	        geometry street_geometry <- i.shape;
                list<point> street_geometry_points_list <- street_geometry.points;
                point first_vertex <- street_geometry.points[0].location;
                point last_vertex <- street_geometry.points[length(street_geometry_points_list) - 1].location;
                graph street_graph <- as_edge_graph(nearest_walkable_street);
     	 
     	//create a not visible walker which creates points along a street (1m steps) - the closest point to the adress will be saved 
     	create walker number: 1 {
            if first_vertex != nil {
                location <- first_vertex; 
            } else {
                write "Error: First vertex is nil!";
                do die;
            }
        }
     		
     	ask walker {
     		bool stillgoing <- true;
     		float min_point_distance <- 1e9;
     		point temp_closest_point_to_vehicle;
     		loop while: stillgoing {
     			do goto target:last_vertex on:street_graph speed: 1.0;
     			float distance_vehicle <- location distance_to myself.location;
     			
     			if (distance_vehicle < min_point_distance) {
            	    min_point_distance <- distance_vehicle;
            	    temp_closest_point_to_vehicle <- self.location;
           	     }
     			if location = last_vertex.location {
     				stillgoing <- false;
     				Possible_points <- Possible_points + temp_closest_point_to_vehicle;
     				}
     			}
     		do die;
     		}     	      	       	     	
     	}   
     	
     	
     	
     	
     	
			//choosing the closest of all 4 closest points
			vehicle_substitute <- with_min_of(Possible_points, each distance_to location);             	
			ask deliverer{
				location <- myself.vehicle_substitute;
				incar <- false;
				total_time <- total_time + enter_exit_time;
				time_for_vehicle_interaction <- time_for_vehicle_interaction+ enter_exit_time;
				parcels <- maxparcels;
				the_target <- nil;
				the_target_name <- nil;
				
			}
		}
	}
	
	//the car drives to the starting position
	action return_to_base {
		
		// find nearest street to set the speed of the vehicle
		float nearest_maxspeed <- 0 #km / #h;
		float min_distance <- 1e9;
		street nearest_street <- nil;
		        
        ask street {
            float distance_to_vehicle <- self.location distance_to (location);
            if (distance_to_vehicle < min_distance) {
                min_distance <- distance_to_vehicle;
                nearest_street <- self;
            }
        }

        if (nearest_street != nil) {
            nearest_maxspeed <- nearest_street.maxspeedT;
        }
		 vspeed <- nearest_maxspeed #km / #h;
		 
		do goto(target: closest_starting_point_on_street, on:the_graph, speed: vspeed);
		time_from_area <- time_from_area + 1;
	ask deliverer {
		location <- myself.location;
		
	}
	}
	   
    

	
	
}
species startpunkt {
	//This species exists only to visualise the starting point of the deliverer and vehicle
	float vspeed;
	rgb color <- #green;
	
	aspect base{
		draw shape color: color border: color;
	} 
}


//creating an experiment
experiment Random_first_parking_spots type: batch repeat:1 until: fertig keep_seed: false{
	list<int> random_numbers <- [];
 init {
 	// it is important to adapt the numbers. loop times will be your sample size and rnd(0,x); x= number of adresses -1
 	loop times:50{
 		random_numbers <- random_numbers + rnd(0, 722);
 		
 	}
 	write "Generated random numbers: " + random_numbers;
 }   

   parameter "Erster Parkplatz" var:simulation_counter among:random_numbers ;

   permanent {
   	
    display Time background: #white refresh: true {
    	chart "Time needed" type: xy {
	        data "Starting point" value: {simulation_counter, total_time};	       
	    }
    }
    display path_length background: #white refresh: true {
    	chart "Length traveled" type: xy {
	        data "Starting point" value: {simulation_counter, total_path_traveled_length};	       
	    }
    }
}

}
