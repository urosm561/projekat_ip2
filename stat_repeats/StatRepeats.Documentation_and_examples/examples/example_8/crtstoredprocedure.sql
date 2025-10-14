CREATE PROCEDURE Insertorupdatefragment(IN textVar varchar(32600), lengthVar INT, out idFrag integer)
LANGUAGE SQL
BEGIN
    IF NOT EXISTS (SELECT ID from Fragment where text = textvar) then
    BEGIN
        INSERT INTO Fragment (text,text_length) VALUES (textVar, lengthVar);
        Values identity_val_local() into idFrag;
    END;
    ELSE SELECT ID INTO idFrag FROM fragment WHERE text = textVar;
    END IF;
 
END 
@
