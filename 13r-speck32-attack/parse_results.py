attack_success_num = 0
for thread_id in range(50):
    with open(f"attack_records_500times/key_recovery_attack_thread_{thread_id}.txt", 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("Kg surviving rate"):
                attack_success_num += 10 * float(line.split()[-1])
                break
print(f"Attack success num is {attack_success_num}.")
print(f"Attack success rate is {attack_success_num / 500}.")