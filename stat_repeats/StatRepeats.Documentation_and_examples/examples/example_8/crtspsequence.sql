CREATE PROCEDURE Sequence(IN namevar varchar(80), versionvar SMALLINT, minlengthvar INT, confidencevar decimal(4,3), typevar varchar(2), probabilityvar INT, out id integer, out donevar integer)
LANGUAGE SQL
BEGIN
DECLARE finvar INT;
set donevar=0;
IF (probabilityvar=1) THEN
       
    	IF NOT EXISTS (SELECT ID from Sequence where name=namevar AND version=versionvar AND min_length=minlengthvar AND confidence=confidencevar AND type=typevar) then
      		BEGIN
        		INSERT INTO Sequence (name,version,min_length,confidence,type,finished) VALUES (namevar,versionvar,minlengthvar,confidencevar,typevar,0);
        		Values identity_val_local() into id;
      		END;
        ELSE 
     		BEGIN
        		SELECT ID,Finished INTO id,finvar FROM Sequence WHERE name=namevar AND version=versionvar AND min_length=minlengthvar AND confidence=confidencevar AND type=typevar;
        		IF (finvar=0) then
          		BEGIN
            		  IF EXISTS (SELECT ID FROM MATCH WHERE seqid=id) then
                               DELETE FROM MATCH WHERE seqid=id;
                          END IF;
          		END; 
        		ELSE set donevar=1;
        		END IF;
   	   	END;    
    	END IF;
ELSE 
        IF NOT EXISTS (SELECT ID from Sequence where name=namevar AND version=versionvar AND min_length=minlengthvar AND confidence IS NULL AND type=typevar) then
      		BEGIN
        	INSERT INTO Sequence (name,version,min_length,type,finished) VALUES (namevar,versionvar,minlengthvar,typevar,0);
        	Values identity_val_local() into id;
      		END;
    	ELSE 
     		BEGIN
        		SELECT ID,Finished INTO id,finvar FROM Sequence WHERE name=namevar AND version=versionvar AND min_length=minlengthvar AND confidence IS NULL AND type=typevar;
        		IF (finvar=0) then
          		BEGIN
            		   IF EXISTS (SELECT ID FROM MATCH WHERE seqid=id) Then
				DELETE FROM MATCH WHERE seqid=id;
                           END IF;
                        END; 
        		ELSE set donevar=1;
        		END IF;
   	   	END;    
    	END IF;


END IF;          
          
 
END 
@