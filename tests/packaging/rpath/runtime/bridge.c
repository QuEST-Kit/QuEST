int external_value(void);
int bridge_value(void) { return external_value() + 1; }
