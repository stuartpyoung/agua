var stageObject = {"project":"Project1","workflow":"abacus","name":"TOPHAT","owner":"admin","number":"1","type":"aligner","username":"syoung","sessionId":"9999999999.9999.999","mode":"removeStage"};
var Workflow = Agua.controllers.workflow.tabPanes[0];

var Stages = Workflow.stages


Agua.removeStage(stageObject);






var stageObject = {"project":"Project1","workflow":"abacus","name":"TOPHAT","owner":"admin","number":"1","type":"aligner","username":"syoung","sessionId":"9999999999.9999.999","mode":"removeStage"};
Agua.removeStage(stageObject);



var _getIndexInArray = function (hasharray, object, keys) {

// GET THE INDEX OF AN OBJECT IN AN ARRAY, IDENTIFY OBJECT USING SPECIFIED KEY VALUES









	if ( hasharray == null )

	{


		return null;

	}



	for ( var i = 0; i < hasharray.length; i++ )

	{

		var arrayObject = hasharray[i];


		var identified = true;

		for ( var j = 0; j < keys.length; j++ )

		{


			if ( arrayObject[keys[j]] != object[keys[j]] )

			{




				identified = false;

				break;

			}

		}




		if ( identified == true )

		{


			return i;

		}

	}



	return null;

};



var _removeObjectFromArray = function (array, object, keys) {

// REMOVE AN OBJECT FROM AN ARRAY, IDENTIFY OBJECT USING SPECIFIED KEY VALUES









	// DEBUG splitace.pl

	var arrayCopy = dojo.clone(array);




	var index = _getIndexInArray(array, object, keys);




	if ( index == null )	return false;




	array.splice(index, 1);




	return true;

};



var _removeArrayFromArray = function(hasharray, removeThis, uniqueKeys) {

// REMOVE AN ARRAY OF OBJECTS FROM A LARGER ARRAY CONTAINING IT,

// IDENTIFYING OBJECTS USING THE SPECIFIED KEY VALUES








	var success = true;

	for ( var j = 0; j < removeThis.length; j++ )

	{






		var removeSuccess = _removeObjectFromArray(hasharray, removeThis[j], uniqueKeys);




		////////// ***** DEBUG ONLY *****

		////////if ( removeSuccess == false )

		////////{				


		////////	success = false;

		////////}

	}




	return success;

};



var stageObject = {"owner":"admin","project":"Project1","workflow":"abacus","appname":"TOPHAT","appnumber":"1"};

var stageParams = dojo.clone(Agua.stageparameters);

var keys = [ "owner", "project", "workflow", "appname", "appnumber" ];
	var values = [ stageObject.owner, stageObject.project, stageObject.workflow, stageObject.appname, stageObject.appnumber ];

var removeTheseStageParams = Agua.filterByKeyValues(dojo.clone(stageParams), keys, values);

var uniqueKeys = [ "username", "project", "workflow", "appname", "appnumber", "name"];

var removeOk = _removeArrayFromArray(stageParams, removeTheseStageParams, uniqueKeys);
